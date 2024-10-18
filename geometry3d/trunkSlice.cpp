
/**
 * @file trunkSlice.cpp
 * @ingroup converters
 *
 * @date 2018/01/11
 *
 *
 *
 * This file is part of the DGtalTools.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <chrono>
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/shapes/MeshVoxelizer.h"
#include "DGtal/kernel/sets/CDigitalSet.h"
#include "DGtal/kernel/domains/HyperRectDomain.h"
#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/io/Display3D.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/images/ImageContainerBySTLVector.h"

#include "CLI11.hpp"

using namespace DGtal;


/**
 @page mesh2vol
 @brief Convert a mesh file into a 26-separated or 6-separated voxelization in a given resolution grid.

@b Usage: mesh2vol [input]

@b Allowed @b options @b are:

@code
positionals:
  1 TEXT:FILE REQUIRED                  mesh file (.off).
  2 TEXT=result.vol                     filename of ouput volumetric file (vol, pgm3d, ...).

Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT:FILE REQUIRED         mesh file (.off).
  -o,--output TEXT=result.vol           filename of ouput volumetric file (vol, pgm3d, ...).
  -m,--margin UINT                      add volume margin around the mesh bounding box.
  -d,--objectDomainBB                   use the digitization space defined from bounding box of input mesh. If seleted, the option --resolution will have no effect.
  -s,--separation UINT:{6,26}=6         voxelization 6-separated or 26-separated.
  -f,--fillValue                                        change the default output  volumetric image value in [1...255].
  -r,--resolution UINT=128              digitization domain size (e.g. 128). The mesh will be scaled such that its bounding box maps to [0,resolution)^3.
@endcode

@b Example:
@code
  $ mesh2vol -i ${DGtal}/examples/samples/tref.off --separation 26 --resolution 256 -o output.vol
@endcode

@see mesh2vol.cpp

*/

using Image2D = ImageContainerBySTLVector < Z2i::Domain, unsigned char >;
using Image3D = ImageContainerBySTLVector < Z3i::Domain, unsigned char >;
using PointR3 = Z3i::RealPoint;
using PointZ3 = Z3i::Point;

Image3D voxelize(const std::string& inputFilename,
                       const std::string& outputFilename,
                       const unsigned int resolution,
                       const unsigned int margin,
                       const unsigned char fillVal)
{
    trace.beginBlock("Preparing the mesh");
    trace.info() << "Reading input file: " << inputFilename;
    Mesh<PointR3> inputMesh;
    inputMesh << inputFilename.c_str();
    trace.info() << " [done]" << std::endl;
    std::pair<PointR3, PointR3> bbox = inputMesh.getBoundingBox();
    trace.info()<< "Mesh bounding box: "<<bbox.first <<" "<<bbox.second<<std::endl;

    const double smax = (bbox.second - bbox.first).max();
    const double factor = resolution / smax;
    const PointR3 translate = -bbox.first;

    trace.info() << "Scale = "<<factor<<" translate = "<<translate<<std::endl;
    
    for(auto it = inputMesh.vertexBegin(), itend = inputMesh.vertexEnd(); it != itend; ++it)
    {
        //scale + translation
        *it += translate;
        *it *= factor;
    }

    //update BB
    bbox = inputMesh.getBoundingBox();

    trace.endBlock();

    trace.beginBlock("Voxelization");
    Z3i::Domain aDomain(bbox.first, bbox.second);
    trace.info()<< "Domain bounding box: "<< aDomain.lowerBound() <<" "<<  aDomain.upperBound() <<std::endl;

    //Digitization step
    Z3i::DigitalSet mySet(aDomain);
    MeshVoxelizer<Z3i::DigitalSet, 6> aVoxelizer;
    aVoxelizer.voxelize(mySet, inputMesh, 1.0);
    trace.info() << " [done] " << std::endl;
    trace.endBlock();

    trace.beginBlock("Exporting");
    // Export the digital set to a vol file
    trace.info()<<aDomain<<std::endl;
    Image3D image(aDomain);

    for(auto p: mySet)
    {
        image.setValue(p, fillVal);
    }
        
    image >> outputFilename.c_str();
    trace.endBlock();

    return image;
}

int main( int argc, char** argv )
{ 
    // parse command line using CLI ----------------------------------------------
    CLI::App app;
    std::string inputFileName;
    std::string outputFileName {"result.vol"};
    unsigned int margin  {0};
    unsigned int separation {6};
    unsigned int resolution {128};
    unsigned char fillValue {128};
    bool unitScale {false};

    app.description("Convert a mesh file into a 26-separated or 6-separated volumetric voxelization in a given resolution grid. \n Example:\n mesh2vol ${DGtal}/examples/samples/tref.off output.vol --separation 26 --resolution 256 ");

    app.add_option("-i,--input,1", inputFileName, "mesh file (.off)." )
        ->required()
        ->check(CLI::ExistingFile);
    app.add_option("-o,--output,2", outputFileName, "filename of ouput volumetric file (vol, pgm3d, ...).",true);
    app.add_option("-m,--margin", margin, "add volume margin around the mesh bounding box.");
    app.add_flag("-d,--objectDomainBB", unitScale, "use the digitization space defined from bounding box of input mesh. If seleted, the option --resolution will have no effect.");
    app.add_option("-f,--fillValue", fillValue, "change the default output  volumetric image value in [1...255].")
        ->expected(0, 255);
    app.add_option("-r,--resolution", resolution,"digitization domain size (e.g. 128). The mesh will be scaled such that its bounding box maps to [0,resolution)^3.", true);


    app.get_formatter()->column_width(40);
    CLI11_PARSE(app, argc, argv);
    // END parse command line using CLI ----------------------------------------------    

    Image3D voxelized_mesh = voxelize(inputFileName, outputFileName, unitScale ? 0 : resolution, margin, fillValue);
  
    return EXIT_SUCCESS;
}

