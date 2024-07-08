/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

/**
 * @file
 * @ingroup geometry3d
 * @author 
 *
 *
 * @date 
 *
 * Source file of the tool
 *
 * This file is part of the DGtal library/DGtalTools-contrib Project.
 */

///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/helpers/Shortcuts.h>
#include <DGtal/helpers/ShortcutsGeometry.h>
#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/io/writers/MeshWriter.h"
#include <DGtal/geometry/tools/RayIntersectionPredicates.h>
#include <DGtal/images/ArrayImageAdapter.h>
#include <DGtal/io/writers/GenericWriter.h>
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/colormaps/BasicColorToScalarFunctors.h"
#include "DGtal/base/FunctorHolder.h"
#include "DGtal/io/colormaps/GradientColorMap.h"

#include <algorithm>
#include <cmath>

#include "CLI11.hpp"

/**

doc :)
**/

using RealPoint     = DGtal::PointVector<3,double>;
using Point         = DGtal::Z2i::Point;
using Domain        = DGtal::Z2i::Domain;
using PolyMesh      = DGtal::PolygonalSurface<RealPoint>;
using VertexRange   = PolyMesh::VertexRange;
using Matrix3       = DGtal::SimpleMatrix<double, 3, 3>;

template<typename TValue>
using Image2D       = DGtal::ImageContainerBySTLVector<Domain, TValue>;


const RealPoint::Component GLOBAL_epsilon = std::numeric_limits<RealPoint::Component>::epsilon();


Matrix3 makeYawRotationMatrix(double aTheta)
{
    double cosTheta = std::cos(aTheta);
    double sinTheta = std::sin(aTheta);
    return Matrix3 {cosTheta,   -sinTheta,  0,
                    sinTheta,   cosTheta,   0,
                    0,          0,          1};
}


struct Ray
{
    // Members
    RealPoint myOrigin;
    RealPoint myDirection;
    
    // Constructor
    Ray()
        : myDirection(1.0, 0.0, 0.0) 
    {}

    Ray(const RealPoint& aOrigin, const RealPoint& aDirection)
        : myOrigin(aOrigin), myDirection(aDirection.getNormalized())
    {}

    // Methods
    double intersectTriangle(const RealPoint& a, const RealPoint& b, const RealPoint& c) const
    {   // Möller–Trumbore algorithm found here : https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
        // could use DGtal::RayIntersectionPredicates instead but it doesn't return the t value (yet) that is interesting to me
        RealPoint v1 = b - a;
        RealPoint v2 = c - a;
        RealPoint ray_cross_v2 = myDirection.crossProduct(v2);
        double det = v1.dot(ray_cross_v2);

        if(std::abs(det) < GLOBAL_epsilon)
        {   // det is too small, ray is considered parallel to the face, we go next
            return -1.0;
        }

        double inv_det = 1.0 / det;
        RealPoint s = myOrigin - a;
        double u = inv_det * s.dot(ray_cross_v2);

        if(u < 0 || u > 1)
        {
            return -1.0;
        }

        RealPoint s_cross_v1 = s.crossProduct(v1);
        double v = inv_det * myDirection.dot(s_cross_v1);

        if(v < 0 || u + v > 1)
        {
            return -1.0;
        }

        double t = inv_det * v2.dot(s_cross_v1);

        if(t > GLOBAL_epsilon)
        {   // intersection exists, replace any previously found intersection if closer to origin
            return t;
        }

        return -1.0;
    }


    std::pair<int, double> intersectSurface(PolyMesh& aPolysurf) const
    {
        PolyMesh::PositionsMap pos = aPolysurf.positions();

        // search for intersect
        size_t i_min = -1;
        double t_min = std::numeric_limits<double>::infinity();
        for(size_t i = 0; i < aPolysurf.nbFaces(); i++)
        {
            VertexRange vertices = aPolysurf.verticesAroundFace(i);

            // check if the number of vertices is acceptable
            if(vertices.size() < 3)
            {   // ignore this face, we need 3 vertices or more
                continue;
            }

            // check for intersetion
            double t = intersectTriangle(pos[vertices[0]], pos[vertices[1]], pos[vertices[2]]);

            if(t > 0 && t < t_min)
            {   // intersection exists, replace any previously found intersection if closer to origin
                t_min = t;
                i_min = i;
            }
        }

        // return face index and dist value
        return std::pair<int, double>(i_min,t_min);
    }
};


struct SampledCenterline
{
    // Members
    std::vector<RealPoint> mySampledPoints;
    double myMinZ, myMaxZ;
    double mySampleSize;

    // Constructors
    SampledCenterline(std::vector<RealPoint> points, int nbSamples)
        : mySampledPoints(nbSamples)
    {
        if(points.size() < 2)
        {
            DGtal::trace.warning() << "Can't have less than 2 points for a centerline." << std::endl;
            return;
        }

        // ensure points are sorted by z coordinate
        std::sort(points.begin(), points.end(), [](RealPoint a, RealPoint b){return a[2] < b[2];});

        myMinZ = points.front()[2];
        myMaxZ = points.back()[2];
        mySampleSize = (myMaxZ - myMinZ) / nbSamples;

        // populate sample list
        size_t j_mem = 0;
        RealPoint p1 = points[0];
        RealPoint p2 = points[0];
        for(size_t i = 0; i < nbSamples; i++)
        {
            double sampleHeight = myMinZ + mySampleSize/2 + mySampleSize * i;

            for(int j = j_mem; j < points.size(); j++)
            {
                p1 = p2;
                p2 = points[j];

                if(p1[2] < sampleHeight
                   && p2[2] >= sampleHeight)
                {   // found the index, interpolating
                    double a = (sampleHeight - p1[2]) / (p2[2] - p1[2]);
                    mySampledPoints[i] = p1 + a * (p2 - p1);

                    // remember j index (no need to go over the already passed points again)
                    j_mem = j;

                    break;
                }
            }
        }
    }

    // Methods
    RealPoint centerlineRepresentant(const RealPoint &p) const 
    {
        unsigned int i = (unsigned int) ceil((p[2]-myMinZ)/mySampledPoints.size());
        assert(i >= 0);
        i = std::min((unsigned int)(mySampledPoints.size()-1), i);
        return mySampledPoints[i];
    }

    RealPoint first()
    {
        return mySampledPoints[0];
    }
};


struct TrunkMapper
{
    // struct to hold the per cell data that we compute
    struct CellData
    {
        // Members
        int myFaceID;
        double myDist;
        RealPoint myNormal;

        // Constructors
        CellData()
            : myFaceID(-1), myDist(0)
        {}

        CellData(int aFaceID, double aDist, const RealPoint& aNormal)
            : myFaceID(aFaceID), myDist(aDist)
        {}
    };

    typedef std::vector< std::vector< CellData > > DataMap;

    // Members
    PolyMesh myTrunkMesh;
    SampledCenterline myTrunkCenter;
    DataMap myDataMap;
    unsigned int myMapWidth;
    unsigned int myMapHeight;

    // Constructors
    TrunkMapper(const PolyMesh& aTrunkMesh, std::vector<RealPoint>& points, unsigned int aWidth, unsigned int aHeight)
        : myDataMap(aHeight, std::vector< CellData >(aWidth)), myTrunkMesh(aTrunkMesh),
          myTrunkCenter(points, aHeight), myMapWidth(aWidth), myMapHeight(aHeight)
    {}

    // Methods
    RealPoint faceBarycenter(int aFaceID)
    {
        VertexRange vertices = myTrunkMesh.verticesAroundFace(aFaceID);

        RealPoint avgPoint;
        for(int vertID : vertices)
        {
            avgPoint += myTrunkMesh.position(vertID);
        }
        return avgPoint/vertices.size();
    }

    double intersectFace(int aFaceID, const Ray& aRay)
    {
        VertexRange vertices = myTrunkMesh.verticesAroundFace(aFaceID);
        RealPoint p1 = myTrunkMesh.position(vertices[0]);
        RealPoint p2 = myTrunkMesh.position(vertices[1]);
        RealPoint p3 = myTrunkMesh.position(vertices[2]);

        // check for intersetion
        return aRay.intersectTriangle(p1, p2, p3);
    }

    CellData navigateMesh(int aFaceID, const Ray& aRay)
    {
        // contains the ID of faces that have been tested for intersection
        std::set<int> visitedFaces;

        // ordered set that contains candidates for intersection, and other data used to compare them
        auto cmp = [](const std::pair<int, double>& a, const std::pair<int, double>& b)
            { return a.second < b.second; };
        std::set< std::pair<int, double>, decltype(cmp)> candidateFaces(cmp);

        bool flag = true;   // flag that we can change to false when we want to stop the search
        while(flag)
        {
            // populate the candidates set with the unvisited neighbours of the last visited face
            for(int vertID : myTrunkMesh.verticesAroundFace(aFaceID))
            {
                for(int faceID : myTrunkMesh.facesAroundVertex(vertID))
                {
                    if(visitedFaces.find(faceID) != visitedFaces.end())
                    {   // face already tested for intersection : skip it
                        continue;
                    }

                    // computing the dot product between the ray's direction and the face's barycenter
                    // the value of the dot product is used to sort the set
                    // higher value -> face is closer (angle wise) to the ray's direction
                    RealPoint faceCenterRN = (faceBarycenter(faceID) - aRay.myOrigin).getNormalized();
                    double dot = faceCenterRN.dot(aRay.myDirection);

                    candidateFaces.insert(std::pair<int, double>(faceID, dot));
                }
            }

            // iterator to the last element of the set (the face with the highest dot product value)
            auto bestCndtIterator = candidateFaces.rbegin();
            if(bestCndtIterator != candidateFaces.rend())
            {   // set has at least an element
                double t = intersectFace(bestCndtIterator->first, aRay);

                if(t > 0)
                {   // intersection, we stop here
                    return CellData(bestCndtIterator->first, t);
                }
                else
                {   // no intersection, we mark the face as visited and do another loop
                    aFaceID = bestCndtIterator->first;                  // last visited face
                    visitedFaces.insert(aFaceID);                       // add it to the visited face list
                    candidateFaces.erase((++bestCndtIterator).base());  // remove it from the candidate faces
                }
            }
            else
            {   // set is empty, no candidates, exiting the loop
                std::cout << "No faces, exiting" << std::endl;
                flag = false;
            }
        }

        // we only ever get here when the search fails
        // so we pay the high price of going through all of the mesh's faces
        std::pair<int, double> res = aRay.intersectSurface(myTrunkMesh);
        return CellData(res.first, res.second,);
    }
    
    void map()
    {
        Matrix3 rotMat(makeYawRotationMatrix(2 * M_PI / myMapWidth));

        // loop over cells
        int previousFaceID = -1;        // holds the previous cell ID, -1 if it's not available
        for(size_t i = 0; i < myTrunkCenter.mySampledPoints.size(); i++)
        {
            Ray ray(myTrunkCenter.mySampledPoints[i], RealPoint(1.0, 0.0, 0.0));

            for(size_t j = 0; j < myMapWidth; j++)
            {
                if(previousFaceID != -1)
                {
                    myDataMap[i][j] = navigateMesh(previousFaceID, ray);
                }
                else
                {
                    std::pair<int, double> res = ray.intersectSurface(myTrunkMesh);
                    
                    myDataMap[i][j].myFaceID = res.first;
                    myDataMap[i][j].myDist = res.second;
                    
                    previousFaceID = res.first;
                }

                previousFaceID = myDataMap[i][j].myFaceID;

                // rotate ray for next point in line
                ray.myDirection = rotMat * ray.myDirection;
            }

            // change previousFaceID to be the direct cell below rather than the last of last line
            previousFaceID = (myDataMap[i][0].myFaceID == -1 ? myDataMap[i][0].myFaceID : previousFaceID);
        }
    }

    void saveDistMap(const std::string& distMapFilename)
    {
        Domain dom(Point(0,0), Point(myMapWidth -1,myMapHeight -1));
        Image2D<double> distMapImage(dom);

        for(int i = 0; i < myMapHeight; i++)
        {
            for(int j = 0; j < myMapWidth; j++)
            {
                double d = myDataMap[i][j].myDist;
                if(std::isfinite(d))
                {
                    distMapImage.setValue(Point(j,i), myDataMap[i][j].myDist);
                }
                else
                {
                    distMapImage.setValue(Point(j,i), 0);
                }
            }
        }

        auto minmax = std::minmax_element(distMapImage.constRange().begin(), distMapImage.constRange().end());

        std::cout << *minmax.first << "\t" << *minmax.second << std::endl;

        DGtal::GradientColorMap<double> distcolormap(*minmax.first, *minmax.second);
        distcolormap.addColor( DGtal::Color( 180, 0, 0 ) );   // blue
	    distcolormap.addColor( DGtal::Color( 200, 10, 10 ) );

        DGtal::STBWriter< Image2D<double>, DGtal::GradientColorMap<double> >::exportPNG(distMapFilename, distMapImage, distcolormap);
    }


    void saveNormalMap(const std::string& normalMapFilename)
    {
        Domain dom(Point(0,0), Point(myMapWidth -1,myMapHeight -1));
        Image2D<DGtal::Color> distMapImage(dom);

        for(int i = 0; i < myMapHeight; i++)
        {
            for(int j = 0; j < myMapWidth; j++)
            {
                double d = myDataMap[i][j].myDist;
                if(std::isfinite(d))
                {
                    // color
                }
                else
                {   // Transparent Color
                    distMapImage.setValue(Point(j,i), DGtal::Color(0,0,0,0));
                }
            }
        }

        auto minmax = std::minmax_element(distMapImage.constRange().begin(), distMapImage.constRange().end());

        std::cout << *minmax.first << "\t" << *minmax.second << std::endl;

        DGtal::STBWriter< Image2D<DGtal::Color> >::exportPNG(normalMapFilename, distMapImage);
    }
};


int main(int argc, char** argv)
{
    std::string meshFilename;
    std::string centerlineFilename;
    std::string outputFilename = "map.png";
    int nbVSamples = 10;
    int nbHSamples = 10;

    DGtal::Mesh<DGtal::Z3i::RealPoint> inputMesh;
    PolyMesh inputPolySurf;

    // parse command line using CLI ----------------------------------------------
    CLI::App app;
    app.description("trunkMeshMap tool to create a map representation of the surface of a tree trunk (mesh).\n");
    
    // inputs
    app.add_option("--inputMesh,1", meshFilename, "an input mesh file in .obj or .off format." )
        ->required()
        ->check(CLI::ExistingFile);
    app.add_option("--inputCenterline,2", centerlineFilename, "an points coordinates input file" )
        ->required()
        ->check(CLI::ExistingFile);
    app.add_option("--nbVertSamples", nbVSamples, "number of vertical samples, also image height in pixels" );
    app.add_option("--nbHoriSamples", nbHSamples, "number of vertical samples, also image height in pixels" );
    
    // outputs
    app.add_option("-o,--output", outputFilename, "an output image file.", true );

    app.get_formatter()->column_width(40);
    CLI11_PARSE(app, argc, argv);

    // END parse command line using CLI ----------------------------------------------

    // read input mesh and transform into PolyMesh
    DGtal::trace.info() << "Reading input mesh...";
    inputMesh << meshFilename;
    DGtal::trace.info() << " [done] (#vertices: " << inputMesh.nbVertex() << ")" << std::endl;
    
    inputMesh.removeIsolatedVertices();

    DGtal::trace.info() << "Mesh into PolyMesh...";
    if( DGtal::MeshHelpers::mesh2PolygonalSurface(inputMesh, inputPolySurf))
    {
        DGtal::trace.info() << " [done] (#vertices: " << inputPolySurf.nbVertices() << ")" << std::endl;
    }
    else
    {
        DGtal::trace.error() << " [failed] (three faces sharing an edge/wrong number of vertices/butterfly neighborhood)" << std::endl;
    }

    DGtal::trace.info() << "Reading input centerline...";
    std::vector<RealPoint> centerline = DGtal::PointListReader<RealPoint>::getPointsFromFile(centerlineFilename);
    DGtal::trace.info() << " [done] (#points: " << centerline.size() << ")" <<  std::endl;

    TrunkMapper TM(inputPolySurf, centerline, nbHSamples, nbVSamples);
    TM.map();

    TM.saveDistMap(outputFilename);
}
