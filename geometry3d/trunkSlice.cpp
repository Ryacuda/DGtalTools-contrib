
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

#include <boost/format.hpp>

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
using SliceImageAdapter = ConstImageAdapter<Image3D,
                                            Image2D::Domain,
                                            functors::Projector< Z3i::Space>,
                                            Image3D::Value,
                                            functors::Identity>;

Image3D voxelizeMesh(Mesh<PointR3>& input_mesh,
                 const unsigned int number_of_slice,
                 const unsigned char fillVal)
{
    
    std::pair<PointR3, PointR3> bbox = input_mesh.getBoundingBox();
    trace.info()<< "Mesh bounding box: "<<bbox.first <<" "<<bbox.second<<std::endl;

    const double smax = (bbox.second - bbox.first).max();
    const double factor = number_of_slice / smax;
    const PointR3 translate = -bbox.first;

    trace.info() << "Scale = "<<factor<<" translate = "<<translate<<std::endl;
    
    for(auto it = input_mesh.vertexBegin(), itend = input_mesh.vertexEnd(); it != itend; ++it)
    {
        //scale + translation
        *it += translate;
        *it *= factor;
    }

    //update BB
    bbox = input_mesh.getBoundingBox();

    trace.endBlock();

    trace.beginBlock("Voxelization");
    Z3i::Domain aDomain(bbox.first, bbox.second);
    trace.info()<< "Domain bounding box: "<< aDomain.lowerBound() <<" "<<  aDomain.upperBound() <<std::endl;

    //Digitization step
    Z3i::DigitalSet mySet(aDomain);
    MeshVoxelizer<Z3i::DigitalSet, 6> aVoxelizer;
    aVoxelizer.voxelize(mySet, input_mesh, 1.0);
    trace.endBlock();

    // Export the digital set to a vol file
    Image3D image(aDomain);

    for(auto p: mySet)
    {
        image.setValue(p, fillVal);
    }

    return image;
}

void sliceVol(const Image3D vol, std::string output_basefilename, unsigned int slice_orientation)
{
    std::string output_fileext = output_basefilename.substr(output_basefilename.find_last_of(".")+1);

    trace.beginBlock("Slicing");

#pragma omp parallel for schedule(dynamic)
    for(size_t i = 0; i < vol.domain().upperBound()[slice_orientation]; i++ )
    {
        trace.info() << "Exporting slice image "<< i ;

        functors::Projector<Z2i::Space>  invFunctor; invFunctor.initRemoveOneDim(slice_orientation);
        Z2i::Domain domain2D(invFunctor(vol.domain().lowerBound()),
                             invFunctor(vol.domain().upperBound()));
        
        functors::Projector<Z3i::Space> aSliceFunctor(i); aSliceFunctor.initAddOneDim(slice_orientation);

        const functors::Identity identityFunctor{};

        SliceImageAdapter sliceImage( vol, domain2D, aSliceFunctor, identityFunctor);

        std::stringstream output_filename;
        output_filename << output_basefilename << "_" <<  boost::format("%|05|")% i <<"."<< output_fileext;

        trace.info() << ": "<< output_filename.str() ;

        GenericWriter<SliceImageAdapter>::exportFile(output_filename.str(), sliceImage);

        trace.info() << " [done]"<< std::endl;
    }

    trace.endBlock();
}

int main( int argc, char** argv )
{ 
    // parse command line using CLI ----------------------------------------------
    CLI::App app;
    std::string input_filename;
    std::string output_basefilename {"result.vol"};
    std::pair<unsigned int, unsigned int> slice_size {192, 192};
    unsigned int separation {6};
    unsigned int number_of_slice {128};
    unsigned char fill_value {128};
    unsigned int slice_axis {2};

    app.description("Convert a mesh file into a 26-separated or 6-separated volumetric voxelization in a given resolution grid. \n Example:\n mesh2vol ${DGtal}/examples/samples/tref.off output.vol --separation 26 --resolution 256 ");

    app.add_option("-i,--input,1", input_filename, "mesh file (.off)." )
        ->required()
        ->check(CLI::ExistingFile);
    app.add_option("-o,--output,2", output_basefilename, "base_name.extension:  extracted 2D slice volumetric files (will result n files base_name_xxx.extension)",true);
    
    app.add_option("-f,--fillValue", fill_value, "change the default output  volumetric image value in [1...255].")
        ->expected(0, 255);
    app.add_option("-s,--sliceSize", slice_size, "desired slice width and height (default 192x192).")
        ->expected(2);
    app.add_option("-a,--sliceAxis", slice_axis, "specify the slice orientation for which the slice are defined (by default =2 (Z direction))", true)
        -> check(CLI::IsMember({0, 1, 2}));
    app.add_option("-n,--nslice", number_of_slice,"digitization domain size (e.g. 128). The mesh will be scaled such that its bounding box largest dim is the number_of_slice", true);


    app.get_formatter()->column_width(40);
    CLI11_PARSE(app, argc, argv);
    // END parse command line using CLI ----------------------------------------------

    trace.beginBlock("Preparing the mesh");
    trace.info() << "Reading input file: " << input_filename;
    Mesh<PointR3> input_mesh;
    input_mesh << input_filename;
    trace.info() << " [done]" << std::endl;

    Image3D voxelized_mesh = voxelizeMesh(input_mesh, number_of_slice, fill_value);

    sliceVol(voxelized_mesh, output_basefilename, slice_axis);

    //voxelized_mesh >> output_filename;
  
    return EXIT_SUCCESS;
}