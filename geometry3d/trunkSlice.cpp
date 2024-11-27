
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
#include <random>

#include "DGtal/helpers/StdDefs.h"
#include "DGtal/shapes/MeshVoxelizer.h"
#include "DGtal/kernel/sets/CDigitalSet.h"
#include "DGtal/kernel/domains/HyperRectDomain.h"
#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/io/Display3D.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/images/ImageContainerBySTLVector.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <boost/algorithm/minmax_element.hpp>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/kernel/sets/DigitalSetInserter.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "DGtal/images/IntervalForegroundPredicate.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"
#include "DGtal/io/colormaps/GrayscaleColorMap.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/io/colormaps/BasicColorToScalarFunctors.h"

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
  -f,--fillValue                        change the default output  volumetric image value in [1...255].
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

template<typename TContainer>
void createDistanceMap(const TContainer& image, const std::string& filename_output, bool invert = true)
{
    Image2D img(image.domain());

    size_t i = 0;
    if(invert)
    {
        for(auto& value : img)
        {
            value = 255 - *(image.constRange().begin() + i);
            i++;
        }
    }

    typedef functors::IntervalForegroundPredicate<Image2D> Binarizer; 
    Binarizer b(img,128, 255); 
    typedef DGtal::DistanceTransformation<Z2i::Space, Binarizer, Z2i::L2Metric> DTL2;
    DTL2 dt(&img.domain(),&b, &Z2i::l2Metric );

    double maxDT = (*boost::first_max_element(dt.constRange().begin(), dt.constRange().end()));

    i = 0;
    for(auto& value : img)
    {
        value = *(dt.constRange().begin() + i);
        i++;
    }

    DGtal::STBWriter<Image2D, DGtal::HueShadeColorMap<unsigned char,2> >::exportPNG(filename_output, img, DGtal::HueShadeColorMap<unsigned char,2>(0, maxDT));
}

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

void sliceVol(const Image3D vol, std::string output_filename, unsigned int slice_orientation)
{
    std::string output_fileext = output_filename.substr(output_filename.find_last_of(".")+1);
    std::string output_basename = output_filename.substr(0, output_filename.find_last_of("."));

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
        std::stringstream dist_filename;
        output_filename << output_basename << "_" <<  boost::format("%|05|")% i <<"."<< output_fileext;
        dist_filename << output_basename << "_dist_" <<  boost::format("%|05|")% i <<".png";

        trace.info() << ": "<< output_filename.str() ;

        GenericWriter<SliceImageAdapter>::exportFile(output_filename.str(), sliceImage);

        createDistanceMap<SliceImageAdapter>(sliceImage, dist_filename.str());

        trace.info() << " [done]"<< std::endl;
    }

    trace.endBlock();
}

class VariationSlicer
{
public:
    // Constructors
    VariationSlicer(const Mesh<PointR3>& baseMesh, const std::vector<PointR3>& pith, bool saveMeshes = false, size_t displacement_nb = 3)
     : myBaseMesh(baseMesh), myPith(pith), saveMesh(saveMeshes)
    {
        if(pith.size() != 0)
        {
            myP1 = PointR3(0, 0, (*pith.begin())[2]);
            myP4 = PointR3(0, 0, (*(pith.end()-1))[2]);
            initBezierPoints(displacement_nb);
        }
    }

    void generateSlices()
    {
        for(size_t i = 0; i < myP2.size(); i++)
        {
            Mesh<PointR3> mesh_copy = myBaseMesh;

            for(auto it = mesh_copy.vertexBegin(), itend = mesh_copy.vertexEnd(); it != itend; ++it)
            {
                PointR3 t = getTranslation(*it, i);
                t[2] = 0;   // don't change z coordinate
                *it += t;
            }

            if(saveMesh)
            {
                mesh_copy >> "mesh_deform_" + std::to_string(i) + ".off";
            }

            Image3D voxelized_mesh = voxelizeMesh(mesh_copy, myNumberOfSlices, myFillValue);

            sliceVol(voxelized_mesh, mySliceBaseFilename, myBaseSliceAxis);
        }
    }

    // Methods
    void initBezierPoints(size_t displacement_nb)
    {
        // generate two more points evenly spaced along z between p1 and p4
        double dz = myP4[2] - myP1[2];

        double amplitude = dz / 2;
        std::random_device rd;  // Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<double> dis(-amplitude, amplitude);

        // create P2 and P3 for each displacement
        for(size_t i = 0; i < displacement_nb; i++)
        {
            myP2.emplace_back(dis(gen), dis(gen), myP1[2] + dz/3);        // p2
            myP3.emplace_back(dis(gen), dis(gen), myP1[2] + 2 * dz/3);    // p3
        }
    }

    void initSliceParameters(unsigned int numberOfSlices, unsigned char fillValue, std::string sliceBaseFilename, unsigned int baseSliceAxis)
    {
        myNumberOfSlices = numberOfSlices;
        myFillValue = fillValue;
        mySliceBaseFilename = sliceBaseFilename;
        myBaseSliceAxis = baseSliceAxis;
    }

    PointR3 getTranslation(const PointR3& p, size_t dspt_index)
    {
        if(dspt_index >= myP2.size())
        {
            std::cout << "displacement index out of bound" << std::endl;
            exit;
        }

        double t = (p[2] - myP1[2]) / (myP4[2] - myP1[2]);
        double oneminust = 1 - t;
        return myP1 * oneminust * oneminust * oneminust +
               myP2[dspt_index] * oneminust * oneminust * t +
               myP3[dspt_index] * oneminust * t * t +
               myP4 * t * t * t;
    }

private:
    // data we work on
    Mesh<PointR3> myBaseMesh;
    std::vector<PointR3> myPith;
    bool saveMesh;

    // deform parameters
    PointR3 myP1;
    std::vector<PointR3> myP2;
    std::vector<PointR3> myP3;
    PointR3 myP4;

    // slice parameters
    unsigned int myNumberOfSlices;
    unsigned char myFillValue;
    std::string mySliceBaseFilename;
    unsigned int myBaseSliceAxis;
};

int main( int argc, char** argv )
{ 
    // parse command line using CLI ----------------------------------------------
    CLI::App app;
    std::string mesh_filename;
    std::string centerline_filename;
    std::string output_basefilename {"result.png"};
    std::pair<unsigned int, unsigned int> slice_size {192, 192};
    unsigned int separation {6};
    unsigned int number_of_slice {128};
    unsigned char fill_value {255};
    unsigned int slice_axis {2};

    app.description("Convert a mesh file into a 26-separated or 6-separated volumetric voxelization in a given resolution grid. \n Example:\n mesh2vol ${DGtal}/examples/samples/tref.off output.vol --separation 26 --resolution 256 ");

    app.add_option("-m,--mesh,1", mesh_filename, "mesh file (.off)." )
        ->required()
        ->check(CLI::ExistingFile);
    app.add_option("-c,--centerline,2", centerline_filename, "mesh file (.off)." )
        ->required()
        ->check(CLI::ExistingFile);
    app.add_option("-o,--output,3", output_basefilename, "base_name.extension:  extracted 2D slice volumetric files (will result n files base_name_xxx.extension)",true);
    
    app.add_option("-f,--fillValue", fill_value, "change the default output  volumetric image value in [1...255].")
        ->expected(0, 255);
    app.add_option("-s,--sliceSize", slice_size, "desired slice width and height (default 192x192).")
        ->expected(2);
    app.add_option("-a,--sliceAxis", slice_axis, "specify the slice orientation for which the slice are defined (by default =2 (Z direction))", true)
        -> check(CLI::IsMember({0, 1, 2}));
    app.add_option("-n,--nslice", number_of_slice,"digitization domain size (e.g. 128). The mesh will be scaled such that its bounding box largest dim is the number_of_slice", true);
    auto save_meshes = app.add_flag("--saveMeshes", "saves the deformed meshes on the disk, disabled by default")
        ->default_val(false);


    app.get_formatter()->column_width(40);
    CLI11_PARSE(app, argc, argv);
    // END parse command line using CLI ----------------------------------------------

    trace.beginBlock("Preparing the mesh");
    trace.info() << "Reading mesh file: " << mesh_filename;
    Mesh<PointR3> input_mesh;
    input_mesh << mesh_filename;
    trace.info() << " [done]" << std::endl;

    trace.info() << "Reading centerline file: " << centerline_filename;
    std::vector<PointR3> pith;
    pith = PointListReader<PointR3>::getPointsFromFile(centerline_filename);
    trace.info() << " [done]" << std::endl;
    trace.endBlock();

    VariationSlicer vs(input_mesh, pith, save_meshes);
    vs.initSliceParameters(number_of_slice, fill_value, output_basefilename, slice_axis);

    vs.generateSlices();
    
    return EXIT_SUCCESS;
}