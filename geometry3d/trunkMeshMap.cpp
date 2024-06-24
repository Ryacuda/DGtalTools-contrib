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

#include <algorithm>

#include "CLI11.hpp"

///////////////////////////////////////////////////////////////////////////////
using namespace DGtal;
///////////////////////////////////////////////////////////////////////////////


/**


**/

typedef PointVector<3,double>               RealPoint;
typedef PolygonalSurface< RealPoint >       PolyMesh;
typedef PolyMesh::VertexRange               VertexRange;

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
    double intersectTriangle(const RealPoint& a, const RealPoint& b, const RealPoint& c)
    {   // Möller–Trumbore algorithm found here : https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
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

RealPoint intersectSurface()
{
    // raycast parameter
    RealPoint origin(0.5,0.5,0);
    RealPoint direction(1,0,0);
    Ray ray(origin, direction);

    // numerical constant
    

    PolyMesh mesh = makeBasicPolyMesh();
    PolyMesh::PositionsMap pos = mesh.positions();

    std::cout << "Number of faces: " << mesh.nbFaces() << std::endl;

    double t_min = std::numeric_limits<double>::infinity();
    size_t face_id = mesh.nbFaces();
    for(size_t i = 0; i < mesh.nbFaces(); i++)
    {
        VertexRange vertices = mesh.verticesAroundFace(i);

        // check if the number of vertices is acceptable
        if(vertices.size() < 3)
        {   // ignore this face, we need 3 vertices or more
            continue;
        }

        // check for intersetion
        double t = ray.intersectTriangle(pos[vertices[0]], pos[vertices[1]], pos[vertices[2]]);

        if(t > 0 && t < t_min)
        {   // intersection exists, replace any previously found intersection if closer to origin
            t_min = t;
            face_id = i;
        }
    }

    if(t_min < std::numeric_limits<double>::infinity())
    {   // if an intersection has been found, return it
        std::cout << face_id << std::endl;
        return ray.myOrigin + ray.myDirection * t_min;
    }
    
    // else return ray origin
    return origin;
}


int main(int argc, char** argv)
{
    std::string inputFilename;
    std::string outputFilename = "map.png";

    DGtal::Mesh<DGtal::Z3i::RealPoint> inputMesh;
    PolyMesh inputPolySurf;

    // parse command line using CLI ----------------------------------------------
    
    // inputs and output
    CLI::App app;
    app.description("trunkMeshMap tool to create a map representation of the surface of a tree trunk (mesh).\n");
    app.add_option("-i,--input,1", inputFilename, "an input mesh file in .obj or .off format." )
    ->required()
    ->check(CLI::ExistingFile);
    app.add_option("-o,--output,2", outputFilename, "an output image file.", true );

    app.get_formatter()->column_width(40);
    CLI11_PARSE(app, argc, argv);

    // END parse command line using CLI ----------------------------------------------

    /* // read input mesh and transform into PolyMesh
    trace.info() << "Reading input mesh...";
    inputMesh << inputFilename;
    trace.info() << " [done] (#vertices: " << inputMesh.nbVertex() << ")" << std::endl;
    
    inputMesh.removeIsolatedVertices();

    trace.info() << "Mesh into PolyMesh...";
    if( DGtal::MeshHelpers::mesh2PolygonalSurface(inputMesh, inputPolySurf))
    {
        trace.info() << " [done] (#vertices: " << inputPolySurf.nbVertices() << ")" << std::endl;
    }
    else
    {
        trace.info() << " [failed] (three faces sharing an edge/wrong number of vertices/butterfly neighborhood)" << std::endl;
    } */

    std::cout << intersectSurface() << std::endl;
}
