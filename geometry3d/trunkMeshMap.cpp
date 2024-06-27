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

#include <algorithm>
#include <cmath>

#include "CLI11.hpp"

/**


**/

typedef DGtal::PointVector<3,double>                RealPoint;
typedef DGtal::PolygonalSurface<RealPoint>          PolyMesh;
typedef PolyMesh::VertexRange                       VertexRange;

const RealPoint::Component GLOBAL_epsilon = std::numeric_limits<RealPoint::Component>::epsilon();

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

        std::cout << "Number of faces: " << aPolysurf.nbFaces() << std::endl;

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
        size_t c = 0;
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
                c++;

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

        std::cout << c << std::endl;
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

PolyMesh makeBasicPolyMesh()
{   // taken from https://www.dgtal.org/doc/stable/moduleHalfEdgeMesh.html#HEM_sec3_1
    PolyMesh mesh;
    mesh.addVertex( RealPoint( 0, 0, 0 ) ); // vertex 0
    mesh.addVertex( RealPoint( 1, 0, 0 ) ); // vertex 1
    mesh.addVertex( RealPoint( 0, 1, 0 ) ); // vertex 2
    mesh.addVertex( RealPoint( 1, 1, 0 ) ); // vertex 3
    mesh.addVertex( RealPoint( 0.5, 0.5, 1 ) ); // vertex 4
    mesh.addTriangle( 0, 1, 4 );            // triangle 0
    mesh.addTriangle( 1, 3, 4 );            // triangle 1
    mesh.addTriangle( 3, 2, 4 );            // triangle 2
    mesh.addTriangle( 2, 0, 4 );            // triangle 3
    mesh.addQuadrangle( 1, 0, 2, 3 );       // quadrangle 4
    bool ok = mesh.build(); // should be true
    
    return mesh;
}


struct TrunkMapper
{
    // Members
    PolyMesh myTrunkMesh;
    SampledCenterline myTrunkCenter;
    unsigned int myMapWidth;
    unsigned int myMapHeight;

    // Constructors
    TrunkMapper(const PolyMesh& aTrunkMesh, std::vector<RealPoint>& points, unsigned int aWidth, unsigned int aHeight)
        : myTrunkMesh(aTrunkMesh), myTrunkCenter(points, aHeight), myMapWidth(aWidth), myMapHeight(aHeight)
    {}

    // Methods
    RealPoint faceBarycenter(int aFaceID)
    {
        VertexRange vertices = myTrunkMesh.verticesAroundFace(aFaceID);

        RealPoint avgPoint;
        for(const RealPoint& p : vertices)
        {
            avgPoint += p;
        }
        return avgPoint/vertices.size();
    }

    double intersectFace(int aFaceID, const Ray& aRay)
    {
        VertexRange vertices = aPolysurf.verticesAroundFace(aFaceID);

        // check for intersetion
        return aRay.intersectTriangle(pos[vertices[0]], pos[vertices[1]], pos[vertices[2]]);
    }

    void navigateMesh(int aFaceID, const RealPoint& aIntersectApprox, const Ray& aRay)
    {
        std::map<int, double> faces;
        faces[aFaceID] = (faceBarycenter(aFaceID) - aIntersectApprox).norm1();

        double hitDist = aRay.intersectTriangle(faces.front().first);
        while(hitDist < 0)
        {
            faces.
        }
    }
    
    void map(const std::string& aOutputFilename)
    {
        Ray ray(myTrunkCenter.first(), RealPoint(1.0, 0.0, 0.0));

        // find the first triangle that the initial ray connects with
        std::pair<int, double> res = ray.intersectSurface(myTrunkMesh);

        if(res.first >= 0)
        {   // ray intersected a face
            std::cout << "hit !" << std::endl;
        }
    }
};


int main(int argc, char** argv)
{
    std::string meshFilename;
    std::string centerlineFilename;
    std::string outputFilename = "map.png";
    int nbVSamples;
    int nbHSamples;

    DGtal::Mesh<DGtal::Z3i::RealPoint> inputMesh;
    PolyMesh inputPolySurf;

    // parse command line using CLI ----------------------------------------------
    CLI::App app;
    app.description("trunkMeshMap tool to create a map representation of the surface of a tree trunk (mesh).\n");
    
    // inputs
    app.add_option("--inputmesh,1", meshFilename, "an input mesh file in .obj or .off format." )
    ->required()
    ->check(CLI::ExistingFile);
    app.add_option("--inputcenterline,2", centerlineFilename, "an points coordinates input file" )
    ->required()
    ->check(CLI::ExistingFile);
    app.add_option("--nbVertSamples", nbVSamples, "number of vertical samples, also image height in pixels" );
    app.add_option("--nbHoriSamples", nbHSamples, "number of vertical samples, also image height in pixels" );
    
    // outputs
    app.add_option("-o,--output", outputFilename, "an output image file.", true );

    app.get_formatter()->column_width(40);
    CLI11_PARSE(app, argc, argv);

    // END parse command line using CLI ----------------------------------------------

    /* // read input mesh and transform into PolyMesh
    DGtal::trace.info() << "Reading input mesh...";
    inputMesh << inputFilename;
    DGtal::trace.info() << " [done] (#vertices: " << inputMesh.nbVertex() << ")" << std::endl;
    
    inputMesh.removeIsolatedVertices();

    DGtal::trace.info() << "Mesh into PolyMesh...";
    if( DGtal::MeshHelpers::mesh2PolygonalSurface(inputMesh, inputPolySurf))
    {
        DGtal::trace.info() << " [done] (#vertices: " << inputPolySurf.nbVertices() << ")" << std::endl;
    }
    else
    {
        DGtal::trace.info() << " [failed] (three faces sharing an edge/wrong number of vertices/butterfly neighborhood)" << std::endl;
    } */

    DGtal::trace.info() << "Reading input pith coordinates...";
    std::vector<RealPoint> centerline = DGtal::PointListReader<RealPoint>::getPointsFromFile(centerlineFilename);
    DGtal::trace.info() << " [done] (#centerline nb of points:" << centerline.size() << ")" <<  std::endl;

    // test mesh (pyramid)
    PolyMesh mesh = makeBasicPolyMesh();

    TrunkMapper(mesh, centerline, nbHSamples, nbVSamples);
}
