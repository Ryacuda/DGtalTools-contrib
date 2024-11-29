/**
 * @file trunkSlice.cpp
 * @ingroup geometry3d
 *
 * @date 2024/11/29
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
 @page trunkSlice
 @brief Deforms and slices a mesh.

@b Allowed @b options @b are:

@code
positionals:
  1 TEXT:FILE REQUIRED                  mesh file (.off).
  2 TEXT:FILE REQUIRED                  centerline file.
  3 TEXT="result.png"                   base_name.extension:  extracted 2D slice volumetric files (will result n files base_name_##_#####.pgm)

Options:
  -h,--help                             Print this help message and exit
  -m,--mesh TEXT:FILE REQUIRED          mesh file (.off).
  -c,--centerline TEXT:FILE REQUIRED    centerline file.
  -o,--output TEXT                      base_name.extension:  extracted 2D slice volumetric files (will result n files base_name_##_#####.pgm)
  -a,--sliceAxis                        specify the slice orientation for which the slice are defined (by default =2 (Z direction))
  -n,--nslice                           number of slices.
  --saveMeshes                          saves the deformed meshes on the disk, disabled by default.
  -v,--variations                       number of different variations (int) and their amplitude (usually between 0 and 1). Default 3 and 0.2.
@endcode

@b Example:
@code
  $ trunkSlice Elm.off Elm_centerline.xyz path/s.pgm -n 40 --saveMeshes -v 4 0.2

@see trunkSlice.cpp

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

/*
 * Create a distance map from an image, and writes it to the disk.
 * Optional parameter to invert the input image before computing distance map.
 */
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

/*
 * Discretizes an input mesh.
 * Resolution is determined by the number of slices wanted out of it.
 */
Image3D voxelizeMesh(Mesh<PointR3>& input_mesh, const unsigned int number_of_slice)
{
    // scale and translate the mesh to fit a bounding box of the desired size
    std::pair<PointR3, PointR3> bbox = input_mesh.getBoundingBox();
    trace.info()<< "Mesh bounding box: "<<bbox.first <<" "<<bbox.second<<std::endl;

    const double smax = (bbox.second - bbox.first).max();
    const double factor = number_of_slice / smax;
    const PointR3 translate = -bbox.first;

    trace.info() << "Scale = "<<factor<<" translate = "<<translate<<std::endl;
    
    for(auto it = input_mesh.vertexBegin(), itend = input_mesh.vertexEnd(); it != itend; ++it)
    {
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

    Image3D image(aDomain);

    for(auto p: mySet)
    {
        image.setValue(p, 255);
    }

    return image;
}

/*
 * Slices a 3D image along an axis.
 * The slices are written on the disk following the base filename parameter.
 */
void sliceVol(const Image3D vol, std::string output_basename, unsigned int slice_orientation)
{
    trace.beginBlock("Slicing");

#pragma omp parallel for schedule(dynamic)
    size_t i = 0;
    for(i = 0; i < vol.domain().upperBound()[slice_orientation]; i++ )
    {
        //trace.info() << "Exporting slice image "<< i ;

        functors::Projector<Z2i::Space>  invFunctor; invFunctor.initRemoveOneDim(slice_orientation);
        Z2i::Domain domain2D(invFunctor(vol.domain().lowerBound()),
                             invFunctor(vol.domain().upperBound()));
        
        functors::Projector<Z3i::Space> aSliceFunctor(i); aSliceFunctor.initAddOneDim(slice_orientation);

        const functors::Identity identityFunctor{};

        SliceImageAdapter sliceImage( vol, domain2D, aSliceFunctor, identityFunctor);

        std::string output_filename = output_basename + "_" + (boost::format("%|05|")% i).str() + ".pgm";
        std::string dist_filename = output_basename + "_dist_" + (boost::format("%|05|")% i).str() + ".png";

        //trace.info() << ": "<< output_filename.str() ;

        GenericWriter<SliceImageAdapter>::exportFile(output_filename, sliceImage);

        createDistanceMap<SliceImageAdapter>(sliceImage, dist_filename);

        //trace.info() << " [done]"<< std::endl;
    }

    trace.info() << i << " slices written : " << output_basename << "_#####.pgm" << std::endl;

    trace.endBlock();
}

/*
 * Class to manage deformation variations, their naming and initiate the slicing
 */
class VariationSlicer
{
public:
    // Constructors
    VariationSlicer(const Mesh<PointR3>& baseMesh,
                    const std::vector<PointR3>& pith,
                    const std::string& baseFilename,
                    size_t displacement_nb = 3,
                    double amplitude_factor = 0.2,
                    bool saveMeshes = false)
     : myBaseMesh(baseMesh), myPith(pith), saveMesh(saveMeshes)
    {
        if(pith.size() != 0)
        {
            myP1 = PointR3(0, 0, (*pith.begin())[2]);
            myP4 = PointR3(0, 0, (*(pith.end()-1))[2]);
            initBezierPoints(displacement_nb, amplitude_factor);
        }

        // strip extension from base file name
        myBaseOutputName = baseFilename.substr(0, baseFilename.find_last_of("."));
    }

    // Methods

    void generateSlices()
    {
        for(size_t i = 0; i < myP2.size(); i++)
        {
            std::string formatted_var_index = (boost::format("%|02|")% i).str();
            Mesh<PointR3> mesh_copy = myBaseMesh;

            for(auto it = mesh_copy.vertexBegin(), itend = mesh_copy.vertexEnd(); it != itend; ++it)
            {
                PointR3 t = getTranslation(*it, i);
                t[2] = 0;   // don't change z coordinate
                *it += t;
            }

            if(saveMesh)
            {
                mesh_copy >> myBaseOutputName + "_mesh_deform_" + formatted_var_index + ".off";
            }

            Image3D voxelized_mesh = voxelizeMesh(mesh_copy, myNumberOfSlices);
            
            std::string base_slice_filename = myBaseOutputName + "_" + formatted_var_index;

            sliceVol(voxelized_mesh, base_slice_filename, myBaseSliceAxis);
        }
    }

    void initBezierPoints(size_t displacement_nb, double amplitude_factor)
    {
        // generate two more points evenly spaced along z between p1 and p4
        double dz = myP4[2] - myP1[2];

        double amplitude = dz * amplitude_factor;
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

    void initSliceParameters(unsigned int numberOfSlices, unsigned int baseSliceAxis)
    {
        myNumberOfSlices = numberOfSlices;
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
               myP2[dspt_index] * 3 * oneminust * oneminust * t +
               myP3[dspt_index] * 3 * oneminust * t * t +
               myP4 * t * t * t;
    }

private:
    // data we work on
    Mesh<PointR3> myBaseMesh;
    std::vector<PointR3> myPith;
    bool saveMesh;
    std::string myBaseOutputName;

    // deform parameters
    PointR3 myP1;
    std::vector<PointR3> myP2;
    std::vector<PointR3> myP3;
    PointR3 myP4;

    // slice parameters
    unsigned int myNumberOfSlices;
    unsigned int myBaseSliceAxis;
};

int main( int argc, char** argv )
{ 
    // parse command line using CLI ----------------------------------------------
    CLI::App app;
    std::string mesh_filename;
    std::string centerline_filename;
    std::string output_basefilename {"result.png"};
    unsigned int separation {6};
    unsigned int number_of_slice {128};
    unsigned int slice_axis {2};
    std::pair<unsigned int, double> variations {3, 0.2};
    bool save_meshes = false;

    app.description("Convert a mesh file into a 26-separated or 6-separated volumetric voxelization in a given resolution grid. \n Example:\n mesh2vol ${DGtal}/examples/samples/tref.off output.vol --separation 26 --resolution 256 ");

    app.add_option("-m,--mesh,1", mesh_filename, "mesh file (.off)." )
        ->required()
        ->check(CLI::ExistingFile);
    app.add_option("-c,--centerline,2", centerline_filename, "centerline file." )
        ->required()
        ->check(CLI::ExistingFile);
    app.add_option("-o,--output,3", output_basefilename, "base_name.extension:  extracted 2D slice volumetric files (will result n files base_name_##_#####.pgm)",true);
    app.add_option("-a,--sliceAxis", slice_axis, "specify the slice orientation for which the slice are defined (by default =2 (Z direction))", true)
        -> check(CLI::IsMember({0, 1, 2}));
    app.add_option("-n,--nslice", number_of_slice,"number of slices.", true);
    app.add_flag("--saveMeshes", save_meshes, "saves the deformed meshes on the disk, disabled by default");
    app.add_option("-v,--variations", variations, "number of different variations (int) and their amplitude (usually between 0 and 1). Default 3 and 0.2.");


    app.get_formatter()->column_width(40);
    CLI11_PARSE(app, argc, argv);
    // END parse command line using CLI ----------------------------------------------

    // open files
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

    // start slicing
    VariationSlicer vs(input_mesh, pith, output_basefilename, variations.first, variations.second, save_meshes);
    vs.initSliceParameters(number_of_slice, slice_axis);

    vs.generateSlices();
    
    return EXIT_SUCCESS;
}