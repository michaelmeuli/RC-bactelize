/*=========================================================================
 *
 *  Copyright MOSAIC Group, ETHZ and MPI-CBG Dresden
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

/*=========================================================================
 * This file is part of the Region Competition segmentation algorithm 
 * described in:
 *
 * Cardinale J, Paul G, Sbalzarini I (2012), "Discrete region com-
 * petition for unknown numbers of connected regions". 
 * Image Processing, IEEE Transactions on, 21(8):3531--3545
 *
 * In order to ensure financial support and allow 
 * further development of this software, please cite above publication in
 * all your documents and manuscripts that made use of this software. 
 * Thanks a lot!
 *
 * Authors:
 * Janick Cardinale
 *=========================================================================*/


#include <stdlib.h>
//#include <time.h>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <string>

#include <boost/program_options/parsers.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/filesystem.hpp>
#include <map>

#include <boost/system/error_code.hpp>

#include "itkImageSource.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkFrontsCompetitionImageFilter.h"

#include "itkCastImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkRegionalMaximaImageFilter.h"
#include "itkRegionalMinimaImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkAndImageFilter.h"
#include "itkRectangularImageSource.h"
#include "itkVector.h"
#include "itkHistogram.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkBlobDetectionImageFilter.h"
#include "itkSphereBitmapImageSource.h"
#include "itkFFTConvolutionImageFilter.h"


#include "itkInvertIntensityImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkFlatStructuringElement.h"

#include "itkStatisticsImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkChangeInformationImageFilter.h"

#include "itkTimeProbe.h"
#include "itkNumericTraits.h"

#include "itkRCPiecewiseConstantGaussianNoiseEnergy.h"
#include "itkRCPiecewiseSmoothGaussianNoiseEnergy.h"
#include "itkRCPiecewiseSmoothSquareDistEnergy.h"
#include "itkRCPiecewiseSmoothLiEnergy.h"
#include "itkRCPiecewiseConstantDeconvolutionSquareDistEnergy.h"
#include "itkRCPiecewiseConstantDeconvolutionGaussianNoiseEnergy.h"
#include "itkRCPiecewiseConstantDeconvolutionPoissonNoiseEnergy.h"
#include "itkRCPiecewiseConstantSquareDistEnergy.h"
#include "itkRCPiecewiseConstantPoissonNoiseEnergy.h"
#include "itkRCPiecewiseSmoothPoissonNoiseEnergy.h"
#include "itkRCBalloonFlow.h"
#include "itkRCVectorFieldFlow.h"

#include "itkRCDiscreteContourLengthApproxEnergy.h"
#include "itkRCKybicKratkyContourLengthApproxEnergy.h"
#include "itkRCConstantOutwardFlow.h"
#include "itkRCShapePrior1DEnergy.h"
#include "itkRCShapePrior23DEnergy.h"
#include "itkRCSpringPotentialEnergy.h"
#include "itkRCExponentialPotentialEnergy.h"
#include "itkRCCenterOfMassDirectedFlow.h"

#include "itkExtractImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkImageIOBase.h"
#include "itkSCIFIOImageIO.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkStatisticsLabelMapFilter.h"
#include "itkStatisticsLabelObject.h"
#include "itkVector.h"

#ifdef USE_GPU
#include "itkRCGPUPiecewiseConstantSquareDistCurvatureRegEnergy.h"
#endif

struct params{
    /// preprocessing parameters
    bool preproc_do_prefilter;
    float preproc_prefilter_sigma;
    std::string preproc_mask_image;
    bool preproc_do_masking;
    float preproc_masking_sigma;
    float preproc_percentile;
    bool preproc_do_normalize;
    bool preproc_discard_lower_percentile;
    bool preproc_discard_upper_percentile;
    float preproc_bgs_scale;
    
    /// initialization parameters
    std::vector<unsigned int> init_rect_boarder;
    float init_max_search_filter_sigma;
    float init_blob_min;
    float init_blob_max;
    unsigned int init_circle_radius;
    std::string init_image;
    std::string init_mode;

    /// segmentation algorithm parameters
    bool seg_allow_fission;
    bool seg_allow_fusion;
    bool seg_allow_handles;
    bool seg_use_fast_evolution;
    bool seg_use_sobolev_gradients;
    unsigned int seg_num_iter;
    bool seg_debug;
    bool verbose;
    
    bool seg_use_gpu;

    /// energy paramteres
    float energy_coeff_data;
    float energy_coeff_length;
    float energy_coeff_balloon;
    float energy_coeff_outward_flow;
    float energy_coeff_shape;
    float energy_coeff_shape_23D;
    float energy_coeff_vec_field;
    float energy_region_merge_ths;
    float energy_particle_interaction_radius;
    float energy_sobolev_gradient_kernel_sigma;
    float energy_local_window_radius;
    float energy_curvature_mask_radius;
    std::vector<float> energy_shape_centroid_coeff;
    
    unsigned int energy_order_of_moments;
    std::string energy_psf_file_name;
    std::string energy_shape_template_file_name;
    std::string energy_vector_field_file_name;
    std::string energy_ext_name;
    std::string energy_int_name;

    /// in/output image parameters
    std::vector<float> image_spacing;
    std::string image_name;
    std::vector<int> channels;


    /// tuning
    float seg_cell_size;
#ifdef USE_GPU
    unsigned int gpu_workgroupo_size;
#endif

    /// operational stuff (not saved in config file)
    std::string parameter_file;
    void setToDefault(const unsigned int, const unsigned int);
    int load(const std::string &filename);
    void save(const std::string &filename);
};

void params::save(const std::string& filename) {
    using boost::property_tree::ptree;
    ptree pt;

    /// preprocessing parameters
    pt.put("params.preproc.preproc_do_prefilter",preproc_do_prefilter);
    pt.put("params.preproc.preproc_prefilter_sigma",preproc_prefilter_sigma);
    pt.put("params.preproc.preproc_do_masking",preproc_do_masking);
    pt.put("params.preproc.preproc_mask_image",preproc_mask_image);
    pt.put("params.preproc.preproc_masking_sigma",preproc_masking_sigma);
    pt.put("params.preproc.preproc_percentile",preproc_percentile);
    pt.put("params.preproc.preproc_do_normalize",preproc_do_normalize);
    pt.put("params.preproc.preproc_discard_lower_percentile",preproc_discard_lower_percentile);
    pt.put("params.preproc.preproc_discard_upper_percentile",preproc_discard_upper_percentile);
    pt.put("params.preproc.preproc_bgs_scale",preproc_bgs_scale);

    
    /// init parametetrs
    for(int vI = 0; vI < init_rect_boarder.size(); vI++) {
        pt.put("params.init.init_rect_boarder." +
            boost::lexical_cast<std::string, int>(vI), init_rect_boarder[vI]);
    }
    pt.put("params.init.init_max_search_filter_sigma",init_max_search_filter_sigma);
    pt.put("params.init.init_blob_min",init_blob_min);
    pt.put("params.init.init_blob_max",init_blob_max);
    pt.put("params.init.init_circle_radius",init_circle_radius);
    pt.put("params.init.init_image",init_image);
    pt.put("params.init.init_mode",init_mode);

    /// segmentation algorithm parameters
    pt.put("params.seg.seg_allow_fission",seg_allow_fission);
    pt.put("params.seg.seg_allow_fusion",seg_allow_fusion);
    pt.put("params.seg.seg_allow_handles",seg_allow_handles);
    pt.put("params.seg.seg_use_fast_evolution",seg_use_fast_evolution);
    pt.put("params.seg.seg_use_sobolev_gradients",seg_use_sobolev_gradients);
    pt.put("params.seg.seg_num_iter",seg_num_iter);
    pt.put("params.seg.seg_use_gpu",seg_use_gpu);
    pt.put("params.seg.seg_debug",seg_debug);
    pt.put("params.seg.verbose", verbose);



    /// energy paramteres
    pt.put("params.energy.energy_coeff_data",energy_coeff_data);
    pt.put("params.energy.energy_coeff_length",energy_coeff_length);
    pt.put("params.energy.energy_coeff_balloon",energy_coeff_balloon);
    pt.put("params.energy.energy_coeff_outward_flow",energy_coeff_outward_flow);
    pt.put("params.energy.energy_coeff_shape",energy_coeff_shape);
    pt.put("params.energy.energy_coeff_shape_23D",energy_coeff_shape_23D);
    pt.put("params.energy.energy_coeff_vec_field", energy_coeff_vec_field);
    pt.put("params.energy.energy_region_merge_ths",energy_region_merge_ths);
    pt.put("params.energy.energy_particle_interaction_radius",energy_particle_interaction_radius);
    pt.put("params.energy.energy_sobolev_gradient_kernel_sigma",energy_sobolev_gradient_kernel_sigma);
    pt.put("params.energy.energy_local_window_radius",energy_local_window_radius);
    pt.put("params.energy.energy_curvature_mask_radius",energy_curvature_mask_radius);
    pt.put("params.energy.energy_psf_file_name",energy_psf_file_name);
    pt.put("params.energy.energy_shape_template_file_name",energy_shape_template_file_name);
    pt.put("params.energy.energy_vector_field_file_name",energy_vector_field_file_name);
    pt.put("params.energy.energy_ext_name",energy_ext_name);
    pt.put("params.energy.energy_int_name",energy_int_name);
    pt.put("params.energy.energy_order_of_moments",energy_order_of_moments);
    for(int vI = 0; vI < energy_shape_centroid_coeff.size(); vI++) {
        pt.put("params.energy.energy_shape_centroid_coeff." +
        boost::lexical_cast<std::string, int>(vI), energy_shape_centroid_coeff[vI]);
    }
                
    /// in/output image parameters
    for(int vI = 0; vI < image_spacing.size(); vI++) {
        pt.put("params.image.image_spacing." +
            boost::lexical_cast<std::string, float>(vI), image_spacing[vI]);
    }
    pt.put("params.image.image_name",image_name);

    pt.put("params.image.channels.bacteria", channels[0]);
    pt.put("params.image.channels.lysosomes", channels[1]);
    pt.put("params.image.channels.cells", channels[2]);

    /// tuning
    pt.put("params.seg.seg_cell_size",seg_cell_size);
#ifdef USE_GPU
    pt.put("params.gpu.gpu_workgroupo_size",gpu_workgroupo_size);
#endif

    boost::property_tree::write_json(filename, pt);
}

int params::load(const std::string& filename){
  boost::property_tree::ptree pt;
  boost::property_tree::read_json(filename, pt);
  try {

    /// preprocessing parameters
    preproc_do_prefilter = pt.get("params.preproc.preproc_do_prefilter",preproc_do_prefilter);
    preproc_prefilter_sigma = pt.get("params.preproc.preproc_prefilter_sigma",preproc_prefilter_sigma);
    preproc_do_masking = pt.get("params.preproc.preproc_do_masking",preproc_do_masking);
    preproc_mask_image = pt.get("params.preproc.preproc_mask_image",preproc_mask_image);
    preproc_masking_sigma = pt.get("params.preproc.preproc_masking_sigma",preproc_masking_sigma);
    preproc_percentile = pt.get("params.preproc.preproc_percentile",preproc_percentile);
    preproc_do_normalize = pt.get("params.preproc.preproc_do_normalize", preproc_do_normalize);
    preproc_discard_lower_percentile = pt.get("params.preproc.preproc_discard_lower_percentile", preproc_discard_lower_percentile);
    preproc_discard_upper_percentile = pt.get("params.preproc.preproc_discard_upper_percentile", preproc_discard_upper_percentile);
    preproc_bgs_scale = pt.get("params.preproc.preproc_bgs_scale", preproc_bgs_scale);

    /// init parametetrs
    for(int vI = 0; vI < init_rect_boarder.size(); vI++) {
      init_rect_boarder[vI] = pt.get("params.init.init_rect_boarder." +
	     boost::lexical_cast<std::string, int>(vI), init_rect_boarder[vI]);
    }
    init_max_search_filter_sigma = pt.get("params.init.init_max_search_filter_sigma",init_max_search_filter_sigma);
    init_blob_min = pt.get("params.init.init_blob_min",init_blob_min);
    init_blob_max = pt.get("params.init.init_blob_max",init_blob_max);
    init_circle_radius = pt.get("params.init.init_circle_radius",init_circle_radius);
    init_image = pt.get("params.init.init_image",init_image);
    init_mode = pt.get("params.init.init_mode",init_mode);

    /// segmentation algorithm parameters
    seg_allow_fission = pt.get("params.seg.seg_allow_fission",seg_allow_fission);
    seg_allow_fusion = pt.get("params.seg.seg_allow_fusion",seg_allow_fusion);
    seg_allow_handles = pt.get("params.seg.seg_allow_handles",seg_allow_handles);
    seg_use_fast_evolution = pt.get("params.seg.seg_use_fast_evolution",seg_use_fast_evolution);
    seg_use_sobolev_gradients = pt.get("params.seg.seg_use_sobolev_gradients",seg_use_sobolev_gradients);
    seg_num_iter = pt.get("params.seg.seg_num_iter",seg_num_iter);
    seg_use_gpu = pt.get("params.seg.seg_use_gpu",seg_use_gpu);
    seg_debug = pt.get("params.seg.seg_debug",seg_debug);
    verbose = pt.get("params.seg.verbose", verbose);


    
    /// energy paramteres
    energy_coeff_data = pt.get("params.energy.energy_coeff_data",energy_coeff_data);
    energy_coeff_length = pt.get("params.energy.energy_coeff_length",energy_coeff_length);
    energy_coeff_balloon = pt.get("params.energy.energy_coeff_balloon",energy_coeff_balloon);
    energy_coeff_outward_flow = pt.get("params.energy.energy_coeff_outward_flow",energy_coeff_outward_flow);
    energy_coeff_shape = pt.get("params.energy.energy_coeff_shape",energy_coeff_shape);
    energy_coeff_shape_23D = pt.get("params.energy.energy_coeff_shape_23D",energy_coeff_shape_23D);
    energy_coeff_vec_field = pt.get("params.energy.energy_coeff_vec_field",energy_coeff_vec_field);
    energy_region_merge_ths = pt.get("params.energy.energy_region_merge_ths",energy_region_merge_ths);
    energy_particle_interaction_radius = pt.get("params.energy.energy_particle_interaction_radius",energy_particle_interaction_radius);
    energy_sobolev_gradient_kernel_sigma = pt.get("params.energy.energy_sobolev_gradient_kernel_sigma",energy_sobolev_gradient_kernel_sigma);
    energy_local_window_radius = pt.get("params.energy.energy_local_window_radius",energy_local_window_radius);
    energy_curvature_mask_radius = pt.get("params.energy.energy_curvature_mask_radius",energy_curvature_mask_radius);
    energy_psf_file_name = pt.get("params.energy.energy_psf_file_name",energy_psf_file_name);
    energy_shape_template_file_name = pt.get("params.energy.energy_shape_template_file_name",energy_shape_template_file_name);
    energy_vector_field_file_name = pt.get("params.energy.energy_vector_field_file_name",energy_vector_field_file_name);
    energy_ext_name = pt.get("params.energy.energy_ext_name",energy_ext_name);
    energy_int_name = pt.get("params.energy.energy_int_name",energy_int_name);
    energy_order_of_moments = pt.get("params.energy.energy_order_of_moments", energy_order_of_moments);
    for(int vI = 0; vI < energy_shape_centroid_coeff.size(); vI++) {
        energy_shape_centroid_coeff[vI] = pt.get("params.energy.energy_shape_centroid_coeff."+
                boost::lexical_cast<std::string,int>(vI),energy_shape_centroid_coeff[vI]);
    }
    /// in/output image parameters
    for(int vI = 0; vI < image_spacing.size(); vI++) {
      image_spacing[vI] = pt.get("params.image.image_spacing." +
	     boost::lexical_cast<std::string, float>(vI), image_spacing[vI]);
    }
    image_name = pt.get("params.image.image_name",image_name);

    channels[0] = pt.get("params.image.channels.bacteria", channels[0]);
    channels[1] = pt.get("params.image.channels.lysosomes", channels[1]);
    channels[2] = pt.get("params.image.channels.cells", channels[2]);

    /// tuning
    seg_cell_size = pt.get("params.seg.seg_cell_size",seg_cell_size);
#ifdef USE_GPU
    gpu_workgroupo_size = pt.get("params.gpu.gpu_workgroupo_size",gpu_workgroupo_size);
#endif


  } catch (boost::property_tree::ptree_error& e) {
    std::cerr << "Error: bad arguement in file: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}

void params::setToDefault(const unsigned int dimension, const unsigned int channelSize) {

    preproc_do_prefilter = false;
    preproc_prefilter_sigma = 1;

    init_mode = "spheres";

    seg_debug = false;
    verbose = false;
    seg_use_gpu = false;
    preproc_do_masking = false;
    preproc_mask_image = "";
    preproc_masking_sigma = 5;
    preproc_do_normalize = true;
    preproc_discard_upper_percentile = false;
    preproc_discard_lower_percentile = false;
    preproc_bgs_scale = 0;
    
    energy_local_window_radius = 7;
    energy_curvature_mask_radius = 4;
    init_max_search_filter_sigma = 5;
    init_blob_max = 50;
    init_blob_min = 10;
    preproc_percentile = 0.8;
    init_circle_radius = 3;
    seg_allow_fission = true;
    seg_allow_fusion = true;
    seg_allow_handles = true;
    seg_use_fast_evolution = false;
    seg_use_sobolev_gradients = false;

    seg_num_iter = 20;
    seg_cell_size = 20;

    
    energy_sobolev_gradient_kernel_sigma = 2;
    energy_coeff_data = 1;
    energy_coeff_length = 0.04;
    energy_coeff_balloon = 0;
    energy_coeff_outward_flow = 0;
    energy_coeff_shape = 0;
    energy_coeff_shape_23D = 0;
    energy_coeff_vec_field = 0;
    energy_region_merge_ths = 0.2;
    energy_particle_interaction_radius = 4;
    energy_ext_name = "pc";
    energy_int_name = "curv";
    energy_order_of_moments = 20;

    for(unsigned int d = 0; d < dimension; d++) {
        energy_shape_centroid_coeff.push_back(0); // only first is 1
        image_spacing.push_back(10);
        init_rect_boarder.push_back(10);
    }
    energy_shape_centroid_coeff[0] = 1;

    for(unsigned int c = 0; c < channelSize; c++) {
        channels.push_back(c);
    }
#ifdef USE_GPU
    gpu_workgroupo_size = 32;
#endif


}

template<typename T>
std::string Real2Str(const T& aN){
    std::stringstream ss;
    ss << aN;
    return ss.str();
}

void setEnergyParamsDescription(
boost::program_options::options_description* aDesc,
        params* aParams) {
    
    namespace po = boost::program_options;
    
    aDesc->add_options()
    ("length,l", po::value<float > (&(aParams->energy_coeff_length))
    ->default_value(aParams->energy_coeff_length,
    Real2Str(aParams->energy_coeff_length)),
    "length term coefficient")
    
    ("balloon,b", po::value<float > (&(aParams->energy_coeff_balloon))
    ->default_value(aParams->energy_coeff_balloon,
    Real2Str(aParams->energy_coeff_balloon)),
    "outward balloon flow coefficient")
    
    ("data", po::value<float > (&(aParams->energy_coeff_data))
    ->default_value(aParams->energy_coeff_data,
    Real2Str(aParams->energy_coeff_data)),
    "data term coefficient")
    
    ("shape", po::value<float > (&(aParams->energy_coeff_shape))
    ->default_value(aParams->energy_coeff_shape,
    Real2Str(aParams->energy_coeff_shape)),
    "shape term coefficient for 1-D moments")
    
    ("shape_23D", po::value<float > (&(aParams->energy_coeff_shape_23D))
    ->default_value(aParams->energy_coeff_shape_23D,
    Real2Str(aParams->energy_coeff_shape_23D)),
    "shape term coefficient for 2/3-D moments")
    
    ("vecfield", po::value<float > (&(aParams->energy_coeff_vec_field))
    ->default_value(aParams->energy_coeff_vec_field,
    Real2Str(aParams->energy_coeff_vec_field)),
    "vector field flow term coefficient")
    
    ("outwardflow,c", po::value<float > (&(aParams->energy_coeff_outward_flow))
    ->default_value(aParams->energy_coeff_outward_flow,
    Real2Str(aParams->energy_coeff_outward_flow)),
    "coefficient of a constant outward flow")
    
    ("merge_ths,t", po::value<float > (&(aParams->energy_region_merge_ths))
    ->default_value(aParams->energy_region_merge_ths,
    Real2Str(aParams->energy_region_merge_ths)),
    "region merging threshold")
    
    ("ir", po::value<float > (&(aParams->energy_particle_interaction_radius))
    ->default_value(aParams->energy_particle_interaction_radius,
    Real2Str(aParams->energy_particle_interaction_radius)),
    "maximal particle interaction radius used for potential based regularization")

    ("cr,R", po::value<float > (&(aParams->energy_curvature_mask_radius))
    ->default_value(aParams->energy_curvature_mask_radius,
    Real2Str(aParams->energy_curvature_mask_radius)),
    "curvature mask radius for curvature based regularizing flows. See Kybic and Kratky.")
    
    ("lr,r",
    po::value<float>(&(aParams->energy_local_window_radius))
    ->default_value(aParams->energy_local_window_radius,
    Real2Str(aParams->energy_local_window_radius)),
    "Radius of mask for local energy support (in world coordinates, see spacing)")
    
    ("energy,e", po::value<std::string > (&(aParams->energy_ext_name))
    ->default_value(aParams->energy_ext_name),
    "external energy type ([pc | pcGauss | pcPoiss | ps | psGauss | psPoisson"
    " | pcDec | pcDecGauss | pcDecPoisson | psLi])")
    
    ("internal,E", po::value<std::string > (&(aParams->energy_int_name))->default_value(aParams->energy_int_name),
    "internal energy type ([curv | manhattenLength | springPotential | expPotential | centerOfMassDirectedFlow])")
    
    ("psf", po::value<std::string > (&(aParams->energy_psf_file_name)), "psf image filename")
    
    ("shape_image", po::value<std::string > (&(aParams->energy_shape_template_file_name)), 
    "template shape image file name")
    
    ("moment_order", po::value<unsigned int> (&(aParams->energy_order_of_moments)), 
    "set the order of moments for shape priors")
    
    ("centroid",po::value< std::vector<float> >(&(aParams->energy_shape_centroid_coeff))->multitoken(),
    "coefficients for the moment information of centroid")
    
    ("vector_field_image", po::value<std::string > (&(aParams->energy_vector_field_file_name)), 
    "vector image file name (vtk-format)");
}

void setPreprocParamsDescription(
boost::program_options::options_description* aDesc,
        params* aParams) {
    
    namespace po = boost::program_options;
    aDesc->add_options()
        ("preproc_sigma", po::value<float>(&(aParams->preproc_prefilter_sigma)),
        //->default_value(aParams->preproc_prefilter_sigma),
        "Sigma (in world coordinates) of Gaussian used to prefilter")

            ("no_normalization", "Do NOT normalize the input image to the [0,1]")

            ("bgs", po::value<float>(&(aParams->preproc_bgs_scale))
            ->default_value(aParams->preproc_bgs_scale),
            "Background subtraction (high-pass filter). Parameter determine spat scale.")

            ("disc_lower", "Saturate/bound the 0.5% lowest intensity pixels.")

            ("disc_upper", "Saturate/bound the 0.5% highest intensity pixels.")

            ("mask_image", po::value<std::string > (&(aParams->preproc_mask_image)),
            "Set a forbidden region using a mask. F = (mask > 0). Giving a filename disables Gaussian masking.")

            ("masking_sigma", po::value<float>(&(aParams->preproc_masking_sigma)),
            //->default_value(aParams->preproc_masking_sigma),
            "Sigma (in world coordinates) of Gaussian used to find a region of interest")

            ("masking_ths", po::value<float>(&(aParams->preproc_percentile)),
            //->default_value(aParams->preproc_percentile),
            "intensity precentile for ROI");
}

void setInitParamsDescription(
boost::program_options::options_description* aDesc,
        params* aParams) {
    
    namespace po = boost::program_options;
    aDesc->add_options()
    ("init_mode", po::value<std::string>(&(aParams->init_mode)),
    "\"spheres\", \"rect\", \"otsu\", \"blob_det\" or \"file\".")
    
    ("init_image", po::value<std::string > (&(aParams->init_image)),
    "init image filename.")
    
    ("init_rb",po::value< std::vector<unsigned int> >(&(aParams->init_rect_boarder))->multitoken(),
    "(vector) Distance of initialization rectangle from image boarder.")
    
    ("init_sr", po::value<unsigned int>(&(aParams->init_circle_radius))
    ->default_value(aParams->init_circle_radius),
    "Sphere radii of initialization spheres.")
    
    ("init_sigma",po::value<float>(&(aParams->init_max_search_filter_sigma))
    ->default_value(aParams->init_max_search_filter_sigma,
    Real2Str(aParams->init_max_search_filter_sigma)),
    "Gaussian sigma of kernel used to find local maxima.")
    
    ("init_blob_min",po::value<float>(&(aParams->init_blob_min))
    ->default_value(aParams->init_blob_min,
    Real2Str(aParams->init_blob_min)),
    "Minimal blob diameter to consider.")
    
    ("init_blob_max",po::value<float>(&(aParams->init_blob_max))
    ->default_value(aParams->init_blob_max,
    Real2Str(aParams->init_blob_max)),
    "Maximal blob diameter to consider.");
}



void setParamsDescription(
        boost::program_options::options_description* aDesc,
        params* aParams, 
        std::string* aOutputfilename,
        std::string* aOutputresultsname) {

    namespace po = boost::program_options;

    aDesc->add_options()
    ("help,h", "Help message")
    
    ("he", "Show energy related parameters and options")
    
    ("hi", "Show initialization related options")
    
    ("hp", "Show preprocessing related options")
    
    ("image", po::value<std::string > (&(aParams->image_name)), "image filename")
    
    ("output_image,o", po::value<std::string>(aOutputfilename),
    "output image file name")

    ("output_results", po::value<std::string>(aOutputresultsname),
    "output results file name")
    
    ("params",po::value<std::string>(&(aParams->parameter_file)),"parameter filename")
    
    ("spacing,s",po::value< std::vector<float> >(&(aParams->image_spacing))->multitoken(),
    //->default_value(vParams.image_spacing),
    "spacing of all input images (overwrites)")

    ("channels",po::value< std::vector<int> >(&(aParams->channels))->multitoken(),
    //->default_value(vParams.channels),
    "assignment of channels: bacteria, lysosomes, cells (default: 0 1 2)")
    
    ("sobolev_sigma",po::value<float >(&(aParams->energy_sobolev_gradient_kernel_sigma)),
    //->default_value(aParams->energy_sobolev_gradient_kernel_sigma),
    "Sobolev gradient kernel sigma (in world coordinates)")
    
    ("no_fusion", //po::value<bool>(&(aParams->seg_allow_fusion)),->default_value(aParams->seg_allow_fusion),
    "Disallow fusion")
    
    ("no_fission",// po::value<bool>(&(aParams->seg_allow_fusion)),->default_value(aParams->seg_allow_fusion),
    "Disallow splits")
    
    ("no_handle",// po::value<bool>(&(aParams->seg_allow_fusion))->default_value(aParams->seg_allow_fusion),
    "Disallow the creation of handles")
    
    ("iter,i", po::value<unsigned int>(&(aParams->seg_num_iter))
    ->default_value(aParams->seg_num_iter),
    "max number of iterations")
    
    ("gpu,g", "Use the GPU for the energy computation.")
    
    ("verbose,v", "Something for the bored user.")
    
    ("debug,d", "debug mode"); 
}


/*
* Dimension
*/
const unsigned int DIMENSION = 3;
const unsigned int CHANNELSIZE = 3;
        
/*
* Typedefs
*/
typedef float RealType;
typedef float InternalPixelType; // has to be real type
typedef unsigned int LabelAbsPixelType; //OR TBinaryPixel
typedef unsigned short OutputPixelType;
typedef int LabelPixelType;

typedef itk::Image<InternalPixelType, 5> InternalImageType5D;
typedef itk::Image<InternalPixelType, 4> InternalImageType4D;
typedef itk::Image<InternalPixelType, DIMENSION> InternalImageType;
typedef itk::Image<OutputPixelType, DIMENSION> OutputImageType;
typedef itk::Image<LabelAbsPixelType, DIMENSION> LabelAbsImageType; //OR TBinaryImage
typedef itk::Image<LabelPixelType, DIMENSION> LabelImageType;
typedef itk::FrontsCompetitionImageFilter
            <InternalImageType, LabelAbsImageType, LabelAbsImageType> SegmenterType;

InternalImageType::Pointer extractchannel(InternalImageType5D::Pointer image5D, int channelnr) {
    typedef itk::ExtractImageFilter<InternalImageType5D, InternalImageType4D> FilterType5;
    FilterType5::Pointer exfilter5 = FilterType5::New();
    exfilter5->InPlaceOn();
    exfilter5->SetDirectionCollapseToSubmatrix();
    InternalImageType5D::RegionType inputRegion5 = image5D->GetLargestPossibleRegion();
    InternalImageType5D::SizeType size5 = inputRegion5.GetSize();
    size5[4] = 0;   
    InternalImageType5D::IndexType start5 = inputRegion5.GetIndex();
    start5[4] = channelnr;  
    InternalImageType5D::RegionType desiredRegion5;
    desiredRegion5.SetSize(size5);
    desiredRegion5.SetIndex(start5);
    exfilter5->SetExtractionRegion(desiredRegion5); 
    exfilter5->SetInput(image5D);

    typedef itk::ExtractImageFilter<InternalImageType4D, InternalImageType> FilterType4;
    FilterType4::Pointer exfilter4 = FilterType4::New();
    exfilter4->InPlaceOn();
    exfilter4->SetDirectionCollapseToSubmatrix();
    exfilter5->Update();
    InternalImageType4D::RegionType inputRegion4 = exfilter5->GetOutput()->GetLargestPossibleRegion();
    InternalImageType4D::SizeType size4 = inputRegion4.GetSize();
    size4[3] = 0;   
    InternalImageType4D::IndexType start4 = inputRegion4.GetIndex();
    start4[3] = 0;  
    InternalImageType4D::RegionType desiredRegion4;
    desiredRegion4.SetSize(size4);
    desiredRegion4.SetIndex(start4);
    exfilter4->SetExtractionRegion(desiredRegion4); 
    exfilter4->SetInput(exfilter5->GetOutput());
    exfilter4->Update();
    return exfilter4->GetOutput();
}


struct objectStruct{
    int label;
    double x;
    double y;
    double z;
    double physicalSize;
    int numberOfPixels;
    double maxDiameter;
    double roundness;
    double mean_lysosome;
    double mean_macrophage;
};




/**
 * The main function. If you think it is a bit too lengthy, you're very welcome
 * to break it down into pieces. 
 * The order is as follows:
 *  * General Typedefs etc
 *  * Parameter/options business
 *  * Preprocessing
 *  * Initialization
 *  * Energies
 *  * Algorithm
 *  * Output
 */
int main(int argc, char** argv) {

    std::cout << "This is region competition v1.0_dev for " << DIMENSION << " dimensions\n";
    std::string vOutputFileName;
    std::string vOutputResultsName;
 
    /* 
     * Program arguments:
     */

    params vParams;
    vParams.setToDefault(DIMENSION, CHANNELSIZE);
    
    namespace po = boost::program_options;
    po::options_description vGeneralDescriptions("General options");
    setParamsDescription(&vGeneralDescriptions, &vParams, &vOutputFileName, &vOutputResultsName);
    po::options_description vEnergyDescriptions("Energy related options");
    setEnergyParamsDescription(&vEnergyDescriptions, &vParams);
    po::options_description vInitDescriptions("Initialization related options");
    setInitParamsDescription(&vInitDescriptions, &vParams);
    po::options_description vPreprocDescriptions("Preprocessing mode related options");
    setPreprocParamsDescription(&vPreprocDescriptions, &vParams);
    
    po::options_description vAllDescriptions("All options");
    vAllDescriptions.add(vGeneralDescriptions)
    .add(vEnergyDescriptions)
    .add(vInitDescriptions)
    .add(vPreprocDescriptions);
    
    po::positional_options_description p_desc;
    p_desc.add("image", 1).add("init_image", 1);//.add("output_image", 1);

    po::variables_map variables_map;

    po::store(po::command_line_parser(argc, argv).options(vAllDescriptions).positional(p_desc).run(), 
	      variables_map);//, po::command_line_style::unix_style ^ po::command_line_style::allow_short);

    po::notify(variables_map);

    if (variables_map.count("help") || argc < 2) {
        std::cout << vGeneralDescriptions << std::endl;
        return 0;
    }

    if (variables_map.count("he") || argc < 2) {
        std::cout << vEnergyDescriptions << std::endl;
        return 0;
    }
    if (variables_map.count("hp") || argc < 2) {
        std::cout << vPreprocDescriptions << std::endl;
        return 0;
    }
    if (variables_map.count("hi") || argc < 2) {
        std::cout << vInitDescriptions << std::endl;
        return 0;
    }
    

    if(variables_map.count("params")) {
        if(vParams.load(vParams.parameter_file) != 0){
            std::cerr << "Reading parameter file failed. Bye bye. " << std::endl;
            return 1;
        }    
        /// values have been overwritten by the parameter file. Now we read them
        /// again to allow for giving single arguments that overwrite arguments
        /// in the config file.
        po::variables_map variables_map;
        po::options_description vDescriptions("Allowed options");
        setParamsDescription(&vDescriptions, &vParams, &vOutputFileName, &vOutputResultsName);
        po::store(po::command_line_parser(argc, argv).options(vDescriptions).positional(p_desc).run(), 
	      variables_map);//, po::command_line_style::unix_style ^ po::command_line_style::allow_short);
        po::notify(variables_map);
    }

    if (variables_map.count("debug")) {
        vParams.seg_debug = true;
    }
    
    if (variables_map.count("gpu") || vParams.seg_use_gpu) {
#ifndef USE_GPU        
        std::cout << "Warning: Not compiled with USE_GPU. Falling back to CPU." << std::endl;
        vParams.seg_use_gpu = false;
#endif
        vParams.seg_use_gpu = true;
    }

    if (variables_map.count("preproc_sigma")) {
        vParams.preproc_do_prefilter = true;
    }
    
    if(variables_map.count("disc_lowest")){
        vParams.preproc_discard_lower_percentile = true;
    }
    
    if(variables_map.count("disc_upper")){
        vParams.preproc_discard_upper_percentile = true;
    }

    if(variables_map.count("mask_image") 
            || variables_map.count("masking_sigma") 
            || variables_map.count("masking_ths")) {
        vParams.preproc_do_masking = true;
    }

    if(variables_map.count("sobolev_sigma")){
        vParams.seg_use_sobolev_gradients = true;
    }

    if(variables_map.count("no_fusion")){
        vParams.seg_allow_fusion = false;
    }

    if(variables_map.count("no_fission")) {
        vParams.seg_allow_fission = false;
    }

    if(variables_map.count("no_handle")) {
        vParams.seg_allow_handles = false;
    }
    
    if(variables_map.count("no_normalization")) {
        vParams.preproc_do_normalize = false;
    }


    // generate filename if not provided by argument 
    if(!variables_map.count("output_image")) { 
        std::stringstream vOutSS;
        vOutSS << "rc_seg_e_" << vParams.energy_ext_name
                << "_ie_" << vParams.energy_int_name
                << "_lc_" << vParams.energy_coeff_length
                << "_cr_" << vParams.energy_curvature_mask_radius
                << "_ir_" << vParams.energy_particle_interaction_radius
                << "_ac_" << vParams.energy_coeff_balloon
                << "_tc_" << vParams.energy_region_merge_ths
                << "_i_" << vParams.seg_num_iter
                << ".tif";
        vOutputFileName = vOutSS.str();
    } 
    

    
    LabelImageType::SizeType vCellSize;

    InternalImageType::SpacingType vSpacingCL;
    vSpacingCL.Fill(1);
    bool vSpacingProvided = false;
    if(variables_map.count("spacing")) {
        vSpacingProvided = true;
        for (unsigned int vD = 0; vD < DIMENSION; vD++) {
            vSpacingCL[vD] = vParams.image_spacing[vD];
            if(itk::NumericTraits<float>::epsilon() >= vSpacingCL[vD]) {
                std::cerr << "Wrong spacing, bye bye." << std::endl;
            }
        }
    }

    for (unsigned int vD = 0; vD < DIMENSION; vD++) {
        vCellSize[vD] = std::max(3, static_cast<int> (ceil(
                vParams.seg_cell_size * vSpacingCL[0] / vSpacingCL[vD])));
    }


    typedef enum InitializationKindType {
        e_fromFile, e_rect, e_spheres, e_otsu, e_blob_det
    } InitializationKindType;
    InitializationKindType vInitKind;


    if (std::string("spheres").compare(vParams.init_mode) == 0) {
        vInitKind = e_spheres;
    } else if (std::string("rect").compare(vParams.init_mode) == 0) {
        vInitKind = e_rect;
    } else if (std::string("file").compare(vParams.init_mode) == 0) {
        vInitKind = e_fromFile;
    } else if (std::string("otsu").compare(vParams.init_mode) == 0) {
        vInitKind = e_otsu;
    } else if (std::string("blob_det").compare(vParams.init_mode) == 0) {
        vInitKind = e_blob_det;
    } else {
        std::cerr << "Wrong initialization mode. Bye bye." << std::endl;
        return 1;
    }


    if(variables_map.count("verbose")){
        vParams.verbose = true;
    }


    /*
     * Start time measurement
     */
    itk::TimeProbe vTimer;
    vTimer.Start();

    /*
     * Read the data image and set image information correctly.
     */
    itk::SCIFIOImageIO::Pointer io = itk::SCIFIOImageIO::New();
    io->DebugOn();
    typedef itk::ImageFileReader<InternalImageType5D> FileReaderType;
    FileReaderType::Pointer vFileReader = FileReaderType::New();
    std::cout << "reader->GetUseStreaming(): " << vFileReader->GetUseStreaming() << std::endl;
    std::cout << "done checking streaming usage" << std::endl;
    vFileReader->SetImageIO(io);
    vFileReader->SetFileName(vParams.image_name);
    typedef itk::StreamingImageFilter<InternalImageType5D, InternalImageType5D> StreamingFilter;
    StreamingFilter::Pointer streamer = StreamingFilter::New();
    streamer->SetInput(vFileReader->GetOutput());
    streamer->SetNumberOfStreamDivisions(4);

    /// Perform an update to be able to read out region-size information etc.
    try {
        vFileReader->Update();
    } catch (itk::ExceptionObject& vE) {
        std::cerr << "Exception caught when reading data image. The error is: " << std::endl;
        std::cerr << vE << std::endl;
    }


    const itk::ImageIOBase * imageIO = vFileReader->GetImageIO();
    itk::ImageIORegion region = imageIO->GetIORegion();
    int regionDim = region.GetImageDimension();
    std::cout << "--== Metadata from ImageIOBase ==--" << std::endl;
    for (int i = 0; i < regionDim; i++) {
      std::cout << "\tDimension " << i + 1 << " Size: "
                << region.GetSize(i) << std::endl;
    }
    for (int i = 0; i < regionDim; i++) {
      if (region.GetSize(i) > 1) {
        std::cout << "\tSpacing " << i + 1 << ": "
                  << imageIO->GetSpacing(i) << std::endl;
      }
    }
    std::cout << "\tByte Order: "
              << imageIO->GetByteOrderAsString(imageIO->GetByteOrder())
              << std::endl;
    std::cout << "\tPixel Stride: " << imageIO->GetPixelStride() << std::endl;
    std::cout << "\tPixel Type: "
              << imageIO->GetPixelTypeAsString(imageIO->GetPixelType())
              << std::endl;
    std::cout << "\tImage Size (in pixels): "
              << imageIO->GetImageSizeInPixels() << std::endl;
    std::cout << "\tPixel Type: "
              << imageIO->GetComponentTypeAsString(imageIO->GetComponentType())
              << std::endl;
    std::cout << "\tRGB Channel Count: "
              << imageIO->GetNumberOfComponents() << std::endl;
    std::cout << "\tNumber of Dimensions: "
              << imageIO->GetNumberOfDimensions() << std::endl;
    
    InternalImageType5D::Pointer image5D =  InternalImageType5D::New();
    image5D = streamer->GetOutput();
    image5D->Update();

    if(vParams.verbose) {
        std::cout << "Bacteria channel: " << vParams.channels[0] << std::endl;
        std::cout << "Lysosome channel: " << vParams.channels[1] << std::endl;
        std::cout << "Phagosome channel: " << vParams.channels[2] << std::endl;
    }

    InternalImageType::Pointer image3Dbacteria = extractchannel(image5D, vParams.channels[0]);

    typedef itk::ChangeInformationImageFilter<InternalImageType> ChangeInformationImageFilterType;
    ChangeInformationImageFilterType::Pointer vChangeDataImgSpacingFilter =
            ChangeInformationImageFilterType::New();
    vChangeDataImgSpacingFilter->SetInput(image3Dbacteria);    
    vChangeDataImgSpacingFilter->SetChangeSpacing(true);
    InternalImageType::SpacingType vSpacingBadHack;              //because --init_mode blob_det doesn't work with --spacing
    vSpacingBadHack.Fill(1);                                     //vSpacingCL replaced with vSpacingBadHack
    vChangeDataImgSpacingFilter->SetOutputSpacing(vSpacingBadHack);

    std::cout << "Original image spacing of image3Dbacteria: "
              << image3Dbacteria->GetSpacing()                     
              << " was overridden to " << vSpacingBadHack << std::endl;
    vChangeDataImgSpacingFilter->Update();
    InternalImageType::Pointer vDataImagePointer = vChangeDataImgSpacingFilter->GetOutput();  
    
    /*
     * PREPROCESSING
     * Discard defect pixel (replace with quantile value)
     */
    if(vParams.preproc_discard_lower_percentile || vParams.preproc_discard_upper_percentile) {
        float v1Quantile, vMedian, v99Quantile;
        typedef itk::Statistics::ScalarImageToHistogramGenerator<InternalImageType>
                HistogramGeneratorType;
        HistogramGeneratorType::Pointer vHistogramGenerator = HistogramGeneratorType::New();
        vHistogramGenerator->SetInput(vDataImagePointer);
        unsigned int vNumberOfBins = 512;
        vHistogramGenerator->SetNumberOfBins(vNumberOfBins);
        vHistogramGenerator->Compute();

        typedef HistogramGeneratorType::HistogramType HistogramType;
        const HistogramType* vHistogram = vHistogramGenerator->GetOutput();
        
        typedef itk::ThresholdImageFilter<InternalImageType> ThresholdImageFilterType;
        typedef ThresholdImageFilterType::Pointer ThresholdImageFilterPtrType;

        if(vParams.preproc_discard_lower_percentile) {
            ThresholdImageFilterPtrType vThresholdFilter = ThresholdImageFilterType::New();
            vThresholdFilter->SetInput(vDataImagePointer);
            float v1Quantile = vHistogram->Quantile(0, 0.002);
            if(vParams.verbose)
                std::cout << "0.1st quant = " << v1Quantile << std::endl;
            vThresholdFilter->SetLower(v1Quantile);
            vThresholdFilter->SetOutsideValue(v1Quantile);            
            vThresholdFilter->Update();
            vDataImagePointer = vThresholdFilter->GetOutput();
        }
        if(vParams.preproc_discard_upper_percentile) {
            ThresholdImageFilterPtrType vThresholdFilter = ThresholdImageFilterType::New();
            vThresholdFilter->SetInput(vDataImagePointer);
            float v99Quantile = vHistogram->Quantile(0, 0.998);
            if(vParams.verbose)
                std::cout << "99.5 quant = " << v99Quantile << std::endl;
            vThresholdFilter->SetUpper(v99Quantile);
            vThresholdFilter->SetOutsideValue(v99Quantile);            
            vThresholdFilter->Update();
            vDataImagePointer = vThresholdFilter->GetOutput();
        }
    }
    
    /*
     * PREPROCESSING
     *  find mask
     */
    typedef itk::SmoothingRecursiveGaussianImageFilter<InternalImageType, InternalImageType>
            RecursiveGaussianFilterType;

    typedef itk::BinaryThresholdImageFilter<InternalImageType, LabelAbsImageType>
            BinaryThresholdFilterType;

    RecursiveGaussianFilterType::Pointer vRecursiveMaskGaussFilter =
      RecursiveGaussianFilterType::New();
    RecursiveGaussianFilterType::SigmaArrayType vMaskSigmaArray;
    BinaryThresholdFilterType::Pointer vBinaryThresholdFilter = BinaryThresholdFilterType::New();

    if (vParams.preproc_do_masking && vParams.preproc_mask_image.compare("") == 0) {
        float vVarMask[DIMENSION];
        for (unsigned int vD = 0; vD < DIMENSION; vD++) {
            vVarMask[vD] = vParams.preproc_masking_sigma * vSpacingCL[0] *
                    vSpacingCL[0] / vSpacingCL[vD]; //std
            vVarMask[vD] *= vVarMask[vD]; //std^2

            vMaskSigmaArray[vD] = vParams.preproc_masking_sigma * vSpacingCL[0] *
                    vSpacingCL[0] / vSpacingCL[vD];
        }

	vRecursiveMaskGaussFilter->SetInput(vDataImagePointer);
	vRecursiveMaskGaussFilter->SetSigmaArray(vMaskSigmaArray);
	vRecursiveMaskGaussFilter->Update();

        typedef itk::Statistics::ScalarImageToHistogramGenerator<InternalImageType>
                HistogramGeneratorType;
        HistogramGeneratorType::Pointer vHistogramGenerator = HistogramGeneratorType::New();
        // vHistogramGenerator->SetInput(vMaskGaussFilter->GetOutput());
        vHistogramGenerator->SetInput(vRecursiveMaskGaussFilter->GetOutput());
//        vHistogramGenerator->SetInput(vFFTMaskGaussFilter->GetOutput());
        unsigned int vNumberOfBins = 255;
        vHistogramGenerator->SetNumberOfBins(vNumberOfBins);
        vHistogramGenerator->Compute();

        typedef HistogramGeneratorType::HistogramType HistogramType;
        const HistogramType* vHistogram = vHistogramGenerator->GetOutput();

        double vQuantile = vHistogram->Quantile(0, vParams.preproc_percentile);
        std::cout << "quantile: " << vQuantile << std::endl;

        vBinaryThresholdFilter->SetLowerThreshold(vQuantile);
        vBinaryThresholdFilter->SetInsideValue(0);
        vBinaryThresholdFilter->SetOutsideValue(1);
        vBinaryThresholdFilter->SetInput(vRecursiveMaskGaussFilter->GetOutput());
    }
    
    LabelAbsImageType::Pointer vForbiddenImagePointer = NULL;
    if(vParams.preproc_do_masking) {
        if(vParams.preproc_mask_image.compare("") != 0) {
            typedef itk::ImageFileReader<LabelAbsImageType> ForbImageFileReaderType;
            ForbImageFileReaderType::Pointer vForbiddenImageReader = 
                    ForbImageFileReaderType::New();
            vForbiddenImageReader->SetFileName(vParams.preproc_mask_image);
            vForbiddenImageReader->Update();
            vForbiddenImagePointer = vForbiddenImageReader->GetOutput();
        } else {
            vForbiddenImagePointer = vBinaryThresholdFilter->GetOutput();
        }
    } 
    
    /*
     * PREPROCESSING
     *      * Background subtraction
     */
    if(vParams.preproc_bgs_scale > 0) {
        std::cout << "Subtracting background...";
        std::flush(std::cout);        
        
        float vRadius = vParams.preproc_bgs_scale /2.0f;
        
        // Get a circular map (TODO, more elegant)
        typedef itk::SphereBitmapImageSource<InternalImageType> SphereImageType;
        SphereImageType::Pointer vSphere = SphereImageType::New();
        float vVolume = 0;
        SphereImageType::RadiusType vSphereRadius;
        
        if(DIMENSION == 2) {
            vVolume = vRadius * vRadius * M_PI;
        } else if(DIMENSION == 3) {
            vVolume = 4.0f / 3.0f * M_PI * vRadius * vRadius * vRadius;
        }
        vSphere->SetBackgroundValue(0);
        vSphere->SetForegroundValue(1.0f/vVolume);
        for(unsigned int vD = 0; vD < DIMENSION; vD++) {
            vSphereRadius[vD] = vRadius; // TODO: scaling
        }
        vSphere->SetRadius(vSphereRadius);
        
        typedef itk::FFTConvolutionImageFilter<InternalImageType, InternalImageType> ConvolutionFilterType;
        ConvolutionFilterType::Pointer vConvolutionFilter = ConvolutionFilterType::New();
        vConvolutionFilter->SetKernelImage(vSphere->GetOutput());
        vConvolutionFilter->SetInput(vDataImagePointer);
        vConvolutionFilter->Update();
        
        typedef itk::StatisticsImageFilter<InternalImageType> StatsFilterType;
        StatsFilterType::Pointer vStats = StatsFilterType::New();
        vStats->SetInput(vConvolutionFilter->GetOutput());
        vStats->Update();
        
        typedef itk::SubtractImageFilter<InternalImageType, InternalImageType> SubtractionFilterType;
        SubtractionFilterType::Pointer vSubtractor = SubtractionFilterType::New();
        vSubtractor->SetInput1(vDataImagePointer);
        vSubtractor->SetInput2(vConvolutionFilter->GetOutput());
        
        typedef itk::AddImageFilter<InternalImageType, InternalImageType> AddImageFilterType;
        AddImageFilterType::Pointer vAdder = AddImageFilterType::New();
        vAdder->SetInput1(vSubtractor->GetOutput());
        vAdder->SetInput2(vStats->GetMean());
        vAdder->Update();
        
        vDataImagePointer = vAdder->GetOutput();
        
        if(vParams.verbose && (vParams.preproc_bgs_scale > 0)) {
           typedef itk::ImageFileWriter<InternalImageType> FileWriterTypeBS;
           FileWriterTypeBS::Pointer vFileWriterBS = FileWriterTypeBS::New();        
           vFileWriterBS->SetFileName("Bgs.nrrd");
           vFileWriterBS->SetInput(vDataImagePointer);       
           try {
               std::cout << "Backgroud substracted image is written to: " << "Bgs.nrrd" << std::endl;
               vFileWriterBS->Update();
           } catch (itk::ExceptionObject & e) {
               std::cerr << "Exception caught after starting pipeline in main():" << std::endl;
               std::cerr << e << std::endl;
               std::cerr << "bye bye." << std::endl;
               return 1;
           }
        }
        
        std::cout << "finished." << std::endl;
    }
    
    /* PREPROCESSING:
     *      * Normalization input image or scale the energy coefficients
     */
    InternalPixelType vMax;
    InternalPixelType vMin;
    if(vParams.preproc_do_masking) {
        typedef itk::LabelStatisticsImageFilter<InternalImageType, LabelAbsImageType>
        LabelStatsImageFilterType;
        LabelStatsImageFilterType::Pointer vLabelStatsFilter = 
        LabelStatsImageFilterType::New();
        vLabelStatsFilter->SetInput(vDataImagePointer);
        vLabelStatsFilter->SetLabelInput(vForbiddenImagePointer);
        vLabelStatsFilter->Update();
        vMax = vLabelStatsFilter->GetMaximum(0);
        vMin = vLabelStatsFilter->GetMinimum(0);
    } else {
        typedef itk::StatisticsImageFilter<InternalImageType> StatisticsImageFilterType;
        StatisticsImageFilterType::Pointer vInpStatsFilter =  StatisticsImageFilterType::New();
        vInpStatsFilter->SetInput(vDataImagePointer);
        vInpStatsFilter->Update();
        vMax = vInpStatsFilter->GetMaximum();
        vMin = vInpStatsFilter->GetMinimum();
    }
    float vScale = (itk::NumericTraits<InternalPixelType>::IsPositive(vMax - vMin)) ? 
        1.0f / (vMax - vMin) : 1.0f / vMax;
    float vInternalEnergyCoeffScaler = (vParams.preproc_do_masking) ? (vMax - vMin) : 1.0f;
    
    typedef itk::ShiftScaleImageFilter<InternalImageType, InternalImageType> ShiftScaleDataFilterType;
    ShiftScaleDataFilterType::Pointer vSSDataFilter = ShiftScaleDataFilterType::New();
    vSSDataFilter->SetInput(vDataImagePointer);
    vSSDataFilter->SetShift(-vMin);
    vSSDataFilter->SetScale(vScale);
    if(vParams.preproc_do_normalize) {
        vDataImagePointer = vSSDataFilter->GetOutput();
    } 
    
    /*
     * PREPROCESSING:
     *  Prefilter the input image (for noisy images)
     */
    typedef itk::DiscreteGaussianImageFilter
    <InternalImageType, InternalImageType> GaussianFilterType; // accurate Gaussian filter
    GaussianFilterType::Pointer vGaussianInputPrefilter = GaussianFilterType::New();
    if (vParams.preproc_do_prefilter) {
        float vVar3[DIMENSION];
        for (unsigned int vD = 0; vD < DIMENSION; vD++) {
            vVar3[vD] = vParams.preproc_prefilter_sigma * vSpacingCL[0] * vSpacingCL[0] / vSpacingCL[vD]; //std
            vVar3[vD] *= vVar3[vD]; //std^2

        }
        vGaussianInputPrefilter->SetVariance(vVar3);
        vGaussianInputPrefilter->SetInput(vDataImagePointer);
    }

    if (vParams.preproc_do_prefilter) {
        vDataImagePointer = vGaussianInputPrefilter->GetOutput();
    } 


    /*
     * INITIALIZATION:
     *  read file or generate rectangle
     */
    typedef itk::ImageFileReader<LabelAbsImageType> InitImageFileReaderType;
    typedef itk::RectangularImageSource<LabelAbsImageType> RectImageSourceType;
    InitImageFileReaderType::Pointer vInitImageFileReader = InitImageFileReaderType::New();
    RectImageSourceType::Pointer vRectImageSource = RectImageSourceType::New();

    if (vInitKind == e_fromFile) {
        vInitImageFileReader->SetFileName(vParams.init_image);
        /// Perform an update to be able to read out region-size information etc.
        try {
            vInitImageFileReader->Update();
        } catch (itk::ExceptionObject & vE) {
            std::cerr << "Exception caught when reading the initialization-image."
                    << " The error is: " << std::endl;
            std::cerr << vE << std::endl;
        }

        // TODO: maybe one should use a changedataimgspacing filter:
        vInitImageFileReader->GetOutput()->SetSpacing(vSpacingCL);
    } else if (vInitKind == e_rect){
        vRectImageSource->SetForegroundValue(1);
        vRectImageSource->SetBackgroundValue(0);

        RectImageSourceType::SizeType vFullImageSize =
                vDataImagePointer->GetBufferedRegion().GetSize();

        RectImageSourceType::OutputImageRegionType vForegroundRegion;
        RectImageSourceType::SizeType vRegionSize;
        RectImageSourceType::IndexType vRegionStartIndex;
        for (unsigned int vD = 0; vD < DIMENSION; vD++) {
            int vS = vFullImageSize[vD] - 2 * vParams.init_rect_boarder[vD];
            if(vS < 1){
                std::cerr << "Initialization rectangle is too large. Bye bye." << std::endl;
                return 1;
            }
            vRegionSize[vD] = vFullImageSize[vD] - 2 * vParams.init_rect_boarder[vD];
            vRegionStartIndex[vD] = vParams.init_rect_boarder[vD];
        }
        vForegroundRegion.SetSize(vRegionSize);
        vForegroundRegion.SetIndex(vRegionStartIndex);
        vRectImageSource->SetSize(vFullImageSize);
        vRectImageSource->SetForegroundRegion(vForegroundRegion);

        vRectImageSource->SetSpacing(vSpacingCL);
    }

    /*
     * INITIALIZATION:
     *  Find initial blobs
     */
    RecursiveGaussianFilterType::Pointer vMaxSearchRecursiveGaussFilter = RecursiveGaussianFilterType::New();
    RecursiveGaussianFilterType::SigmaArrayType vMaxSearchSigmaArray;

    float vVar2[DIMENSION];
    for (unsigned int vD = 0; vD < DIMENSION; vD++) {
        vVar2[vD] = vParams.init_max_search_filter_sigma * vSpacingCL[0] * vSpacingCL[0] / vSpacingCL[vD];//std
        vVar2[vD] *= vVar2[vD]; //std^2

        vMaxSearchSigmaArray[vD] =  vParams.init_max_search_filter_sigma * vSpacingCL[0] * vSpacingCL[0] / vSpacingCL[vD];
    }
    vMaxSearchRecursiveGaussFilter->SetSigmaArray(vMaxSearchSigmaArray);
    vMaxSearchRecursiveGaussFilter->SetInput(vDataImagePointer);

    typedef itk::RegionalMaximaImageFilter<InternalImageType, LabelAbsImageType>
            RegionalMaxImageFilterType;
    RegionalMaxImageFilterType::Pointer vMaxFilter = RegionalMaxImageFilterType::New();
    vMaxFilter->SetInput(vMaxSearchRecursiveGaussFilter->GetOutput());
    vMaxFilter->SetBackgroundValue(0);
    vMaxFilter->SetForegroundValue(1);
    vMaxFilter->SetFullyConnected(false);
    typedef itk::BinaryBallStructuringElement<LabelAbsPixelType, DIMENSION>
            StructuringElementType;
    typedef itk::BinaryDilateImageFilter<LabelAbsImageType, LabelAbsImageType,
            itk::FlatStructuringElement<DIMENSION> >
            DilateImageFilterType;
    DilateImageFilterType::Pointer vDilateFilter = DilateImageFilterType::New();
    vDilateFilter->SetInput(vMaxFilter->GetOutput());
    vDilateFilter->SetDilateValue(1);
    StructuringElementType vStructuringElement;
    vStructuringElement.SetRadius(vParams.init_circle_radius);
    vStructuringElement.CreateStructuringElement();
    itk::FlatStructuringElement<DIMENSION>::RadiusType vSERadius;
    for (unsigned int vD = 0; vD < DIMENSION; vD++) {
        vSERadius[vD] = vParams.init_circle_radius;
    }
    vDilateFilter->SetKernel(itk::FlatStructuringElement<DIMENSION>::Ball(vSERadius));



    /* 
     * INITIALIZATION
     *  Set the init-image pointer according to the users choice
     */
    LabelAbsImageType::Pointer vInitImagePointer;
    if(vInitKind == e_rect) {
        vInitImagePointer = vRectImageSource->GetOutput();
    } else if (vInitKind == e_fromFile) {
        vInitImagePointer = vInitImageFileReader->GetOutput();
    } else if(vInitKind == e_spheres) {
        vDilateFilter->Update(); // ?
        vInitImagePointer = vDilateFilter->GetOutput();
    } else if(vInitKind == e_otsu) {
        typedef itk::RegionalMinimaImageFilter<InternalImageType, LabelAbsImageType>
            RegionalMinImageFilterType;
        RegionalMinImageFilterType::Pointer vMinFilter = RegionalMinImageFilterType::New();
        vMinFilter->SetBackgroundValue(1);
        vMinFilter->SetForegroundValue(0);
        vMinFilter->SetInput(vDataImagePointer);
        vMinFilter->SetFullyConnected(false);
    
        typedef itk::OtsuThresholdImageFilter<InternalImageType, LabelAbsImageType> 
        OtsuThsImageFilterType;
        OtsuThsImageFilterType::Pointer vOtsuThsFilter = OtsuThsImageFilterType::New();
        vOtsuThsFilter->SetInput(vDataImagePointer);
        vOtsuThsFilter->SetInsideValue(0);
        vOtsuThsFilter->SetOutsideValue(1);
        
        typedef itk::AndImageFilter<LabelAbsImageType, LabelAbsImageType, LabelAbsImageType> 
        AndImageFilterType;
        AndImageFilterType::Pointer vAndFilter = AndImageFilterType::New();
        vAndFilter->SetInput1(vMinFilter->GetOutput());
        vAndFilter->SetInput2(vOtsuThsFilter->GetOutput());
        vAndFilter->Update();
        vInitImagePointer = vAndFilter->GetOutput();
    } else if(vInitKind == e_blob_det) {
        typedef itk::BlobDetectionImageFilter<InternalImageType, LabelAbsImageType>
        BlobDetectionFilterType;
        BlobDetectionFilterType::Pointer vBlobDetFilter = BlobDetectionFilterType::New();
        vBlobDetFilter->SetNumberOfScales(5);
        vBlobDetFilter->SetMinDiameter(vParams.init_blob_min);
        vBlobDetFilter->SetMaxDiameter(vParams.init_blob_max);
        vBlobDetFilter->SetInput(vDataImagePointer);
        vBlobDetFilter->Update();
        vInitImagePointer = vBlobDetFilter->GetOutput();
    }


    // Michael Meuli: only to find optimal parameters
    if(vParams.verbose && (vInitKind == e_spheres)) {
       typedef itk::ImageFileWriter<InternalImageType> GaussFileWriterType;
       GaussFileWriterType::Pointer vGaussFileWriter = GaussFileWriterType::New();
       vGaussFileWriter->SetFileName("Gauss.nrrd");
       vGaussFileWriter->SetInput(vMaxSearchRecursiveGaussFilter->GetOutput());
       try {
           std::cout << "Output image is written to: " << "Gauss.nrrd" << std::endl;
           vGaussFileWriter->Update();
       } catch (itk::ExceptionObject & e) {
           std::cerr << "Exception caught after starting pipeline in main():" << std::endl;
           std::cerr << e << std::endl;
           std::cerr << "bye bye." << std::endl;
           return 1;
       }
     }

     if(vParams.verbose && (vInitKind == e_spheres)) {
        typedef itk::ImageFileWriter<LabelAbsImageType> FileWriterTypeSpheres;
        FileWriterTypeSpheres::Pointer vFileWriterSpheres = FileWriterTypeSpheres::New();
        vFileWriterSpheres->SetFileName("Spheres.nrrd");     
        vFileWriterSpheres->SetInput(vDilateFilter->GetOutput());
        try {
            std::cout << "Output image is written to: " << "Spheres.nrrd" << std::endl;
            vFileWriterSpheres->Update();
        } catch (itk::ExceptionObject & e) {
            std::cerr << "Exception caught after starting pipeline in main():" << std::endl;
            std::cerr << e << std::endl;
            std::cerr << "bye bye." << std::endl;
            return 1;
        }
     }

     if(vParams.verbose && (vInitKind == e_blob_det)) {
        typedef itk::ImageFileWriter<LabelAbsImageType> FileWriterTypeBlob;
        FileWriterTypeBlob::Pointer vFileWriterBlob = FileWriterTypeBlob::New();        
        vFileWriterBlob->SetFileName("Blobs.nrrd");
        vFileWriterBlob->SetInput(vInitImagePointer);       
        try {
            std::cout << "Blobs image is written to: " << "Blobs.nrrd" << std::endl;
            vFileWriterBlob->Update();
        } catch (itk::ExceptionObject & e) {
            std::cerr << "Exception caught after starting pipeline in main():" << std::endl;
            std::cerr << e << std::endl;
            std::cerr << "bye bye." << std::endl;
            return 1;
        }
     }
    


    /*
     * ENERGIES:
     *      * All the shape template stuff
     */
    typedef LabelAbsImageType ShapeReferenceImageType;
    typedef itk::ImageFileReader<ShapeReferenceImageType> ShapeTemplateFileReaderType;
    ShapeTemplateFileReaderType::Pointer vShapeTemplateReader =
            ShapeTemplateFileReaderType::New();
    vShapeTemplateReader->SetFileName(vParams.energy_shape_template_file_name);
    if (vParams.energy_coeff_shape > 0|| vParams.energy_coeff_shape_23D > 0 ) {
        try {
            vShapeTemplateReader->Update();
        } catch (itk::ExceptionObject & vE) {
            std::cerr << "Exception caught when reading the shape-template image."
                    << " The error is: " << std::endl;
            std::cerr << vE << std::endl;
        }
    }
    
    /*
     * ENERGIES:
     *      * Vector field image reading
     */
    typedef itk::Vector<float, DIMENSION> VectorType;
    typedef itk::Image<VectorType, DIMENSION> VectorImageType;
    typedef itk::ImageFileReader<VectorImageType> VectorFileReaderType;
    VectorFileReaderType::Pointer vVecFileReader = VectorFileReaderType::New();
    vVecFileReader->SetFileName(vParams.energy_vector_field_file_name);
    if(vParams.energy_coeff_vec_field != 0) {
        try {
            vVecFileReader->Update();
        } catch (itk::ExceptionObject & vE) {
            std::cerr << "Exception caught when reading the vector field image."
                    << " The error is: " << std::endl;
            std::cerr << vE << std::endl;
        }
    }
    
    /*
     * ENERGIES:
     *  All the PSF stuff: Read the file, calc statistics, normalize it
     */
    typedef itk::ImageFileReader<InternalImageType> PSFFileReaderType;
    PSFFileReaderType::Pointer vPSFFileReader = PSFFileReaderType::New();
    vPSFFileReader->SetFileName(vParams.energy_psf_file_name);
    typedef itk::StatisticsImageFilter<InternalImageType> StatisticsFilterType;
    if (variables_map.count("psf")){
        vPSFFileReader->Update();
    }
    vPSFFileReader->GetOutput()->SetSpacing(vSpacingCL);

    /*
     * ENERGIES:
     *      * generate a PSF (TODO)
     */
//                if (m_UseGaussianPSF) {
//                // Here, no PSF has been set by the user. Hence, a Gaussian
//                // approximation is used.
//                typename GaussianImageSourceType::ArrayType vMean;
//                typename GaussianImageSourceType::ArrayType vSigma;
//                typename GaussianImageSourceType::SizeType vGaussImageSize;
//                for (unsigned int vD = 0; vD < m_Dim; vD++) {
//                    vSigma[vD] = m_SigmaPSF[vD];
//                    vGaussImageSize.Fill(m_SigmaPSF[vD] * 6 + 1);
//                }
//                m_GaussianImageSource = GaussianImageSourceType::New();
//                m_GaussianImageSource->SetMean(vMean);
//                m_GaussianImageSource->SetSigma(vSigma);
//                m_GaussianImageSource->SetSize(vGaussImageSize);
//                m_GaussianImageSource->SetNormalized(true);
//                m_GaussianImageSource->Update();
//                m_PSF = m_GaussianImageSource->GetOutput();
//                WriteInputImageTypeToFile("generatedPSF.tif", m_PSF, 100000);
//            }
    
    /* 
     * ENERGIES:
     *      *  Setting up all energies available:
     */
    typedef itk::RCPiecewiseConstantSquareDistEnergy<LabelImageType, InternalImageType> PCSqDistExternalEnergyType;
    PCSqDistExternalEnergyType::Pointer vPCSqDistEnergy = PCSqDistExternalEnergyType::New();
    vPCSqDistEnergy->SetCoefficient(vParams.energy_coeff_data);
    vPCSqDistEnergy->SetRegionMergingThreshold(vParams.energy_region_merge_ths);
    
    typedef itk::RCPiecewiseConstantGaussianNoiseEnergy<LabelImageType, InternalImageType> PCGaussExternalEnergyType;
    PCGaussExternalEnergyType::Pointer vPCGaussEnergy = PCGaussExternalEnergyType::New();
    vPCGaussEnergy->SetCoefficient(vParams.energy_coeff_data);
    vPCGaussEnergy->SetRegionMergingThreshold(vParams.energy_region_merge_ths);
    
    typedef itk::RCPiecewiseConstantPoissonNoiseEnergy<LabelImageType, InternalImageType> PCPoissonExternalEnergyType;
    PCPoissonExternalEnergyType::Pointer vPCPoissonEnergy = PCPoissonExternalEnergyType::New();
    vPCPoissonEnergy->SetCoefficient(vParams.energy_coeff_data);
    vPCPoissonEnergy->SetRegionMergingThreshold(vParams.energy_region_merge_ths);
    
    typedef itk::RCPiecewiseConstantDeconvolutionSquareDistEnergy<LabelImageType, InternalImageType> 
    PCSquaredDistDeconvExternalEnergyType;
    PCSquaredDistDeconvExternalEnergyType::Pointer vPCDecSquaredDistEnergy = PCSquaredDistDeconvExternalEnergyType::New();
    vPCDecSquaredDistEnergy->SetCoefficient(vParams.energy_coeff_data);
    vPCDecSquaredDistEnergy->SetRegionMergingThreshold(vParams.energy_region_merge_ths);
    vPCDecSquaredDistEnergy->SetPSF(vPSFFileReader->GetOutput());
    vPCDecSquaredDistEnergy->SetDoNormalizePSF(true);
//    if(!vParams.use_mcmc) {
    vPCDecSquaredDistEnergy->SetOptimizationMode(true);
//    }
    
    typedef itk::RCPiecewiseConstantDeconvolutionPoissonNoiseEnergy<LabelImageType, InternalImageType> 
    PCPoissonDeconvExternalEnergyType;
    PCPoissonDeconvExternalEnergyType::Pointer vPCDecPoissonEnergy = PCPoissonDeconvExternalEnergyType::New();
    vPCDecPoissonEnergy->SetCoefficient(vParams.energy_coeff_data);
    vPCDecPoissonEnergy->SetRegionMergingThreshold(vParams.energy_region_merge_ths);
    vPCDecPoissonEnergy->SetPSF(vPSFFileReader->GetOutput());
    vPCDecPoissonEnergy->SetDoNormalizePSF(true);
    //    if(!vParams.use_mcmc) {
    vPCDecPoissonEnergy->SetOptimizationMode(true);
//    }
    
    typedef itk::RCPiecewiseConstantDeconvolutionGaussianNoiseEnergy<LabelImageType, InternalImageType> 
    PCGaussDeconvExternalEnergyType;
    PCGaussDeconvExternalEnergyType::Pointer vPCDecGaussEnergy = PCGaussDeconvExternalEnergyType::New();
    vPCDecGaussEnergy->SetCoefficient(vParams.energy_coeff_data);
    vPCDecGaussEnergy->SetRegionMergingThreshold(vParams.energy_region_merge_ths);
    vPCDecGaussEnergy->SetPSF(vPSFFileReader->GetOutput());
    vPCDecGaussEnergy->SetDoNormalizePSF(true);
    //    if(!vParams.use_mcmc) {
    vPCDecGaussEnergy->SetOptimizationMode(true);
//    }
    
    typedef itk::RCPiecewiseSmoothGaussianNoiseEnergy<LabelImageType, InternalImageType> PSGaussExternalEnergyType;
    PSGaussExternalEnergyType::Pointer vPSGaussEnergy = PSGaussExternalEnergyType::New();
    vPSGaussEnergy->SetRadius(vParams.energy_local_window_radius);
    vPSGaussEnergy->SetRegionMergingThreshold(vParams.energy_region_merge_ths);
    vPSGaussEnergy->SetCoefficient(vParams.energy_coeff_data);

    typedef itk::RCPiecewiseSmoothPoissonNoiseEnergy<LabelImageType, InternalImageType> PSPoissonExternalEnergyType;
    PSPoissonExternalEnergyType::Pointer vPSPoissonEnergy = PSPoissonExternalEnergyType::New();
    vPSPoissonEnergy->SetRadius(vParams.energy_local_window_radius);
    vPSPoissonEnergy->SetRegionMergingThreshold(vParams.energy_region_merge_ths);
    vPSPoissonEnergy->SetCoefficient(vParams.energy_coeff_data);
    
    typedef itk::RCPiecewiseSmoothSquareDistEnergy<LabelImageType, InternalImageType> PSSqDistExternalEnergyType;
    PSSqDistExternalEnergyType::Pointer vPSSqDistEnergy = PSSqDistExternalEnergyType::New();
    vPSSqDistEnergy->SetRadius(vParams.energy_local_window_radius);
    vPSSqDistEnergy->SetRegionMergingThreshold(vParams.energy_region_merge_ths);
    vPSSqDistEnergy->SetCoefficient(vParams.energy_coeff_data);
    
    typedef itk::RCPiecewiseSmoothLiEnergy<LabelImageType, InternalImageType> PSLiExternalEnergyType;
    PSLiExternalEnergyType::Pointer vPSLiEnergy = PSLiExternalEnergyType::New();
    vPSLiEnergy->SetCoefficient(vParams.energy_coeff_data);
    vPSLiEnergy->SetRegionMergingThreshold(vParams.energy_region_merge_ths);
    vPSLiEnergy->SetSigma(vParams.energy_local_window_radius);

    typedef itk::RCBalloonFlow<LabelImageType, InternalImageType> BalloonFlowType;
    BalloonFlowType::Pointer vBalloonFlow = BalloonFlowType::New();
    vBalloonFlow->SetCoefficient(vParams.energy_coeff_balloon);
    
    typedef itk::RCVectorFieldFlow<LabelImageType, InternalImageType, VectorImageType> VectorFieldFlowType;
    VectorFieldFlowType::Pointer vVectorFieldFlow = VectorFieldFlowType::New();
    vVectorFieldFlow->SetCoefficient(vParams.energy_coeff_vec_field);
    vVectorFieldFlow->SetVectorImage(vVecFileReader->GetOutput());
    
    typedef itk::RCDiscreteContourLengthApproxEnergyEnergy<LabelImageType> DiscrLengthApproxEnergyType;
    DiscrLengthApproxEnergyType::Pointer vDiscrLengthEnergy = DiscrLengthApproxEnergyType::New();
    vDiscrLengthEnergy->SetCoefficient(vInternalEnergyCoeffScaler * vParams.energy_coeff_length);
    
    typedef itk::RCConstantOutwardFlowEnergy<LabelImageType> ConstOutwardFlowType;
    ConstOutwardFlowType::Pointer vConstOutwardFlow = ConstOutwardFlowType::New();
    vConstOutwardFlow->SetCoefficient(vInternalEnergyCoeffScaler * vParams.energy_coeff_outward_flow);
    
    typedef itk::RCKybicKratkyContourLengthApproxEnergyEnergy<LabelImageType> KKLengthApproxEnergyType;
    KKLengthApproxEnergyType::Pointer vKKLengthEnergy = KKLengthApproxEnergyType::New();
    vKKLengthEnergy->SetCoefficient(vInternalEnergyCoeffScaler * vParams.energy_coeff_length);
    
    typedef itk::RCSpringPotentialEnergy<LabelImageType, SegmenterType> SpringPotEnergyType;
    SpringPotEnergyType::Pointer vSpringPotEnergy = SpringPotEnergyType::New();
    vSpringPotEnergy->SetCoefficient(vInternalEnergyCoeffScaler * vParams.energy_coeff_length); // TODO?
    vSpringPotEnergy->SetParticleInteractionRadius(vParams.energy_particle_interaction_radius);
    
    typedef itk::RCExponentialPotentialEnergy<LabelImageType, SegmenterType> ExponentialPotEnergyType;
    ExponentialPotEnergyType::Pointer vExponentialPotEnergy = ExponentialPotEnergyType::New();
    vExponentialPotEnergy->SetCoefficient(vInternalEnergyCoeffScaler * vParams.energy_coeff_length); // TODO?
    vExponentialPotEnergy->SetParticleInteractionRadius(vParams.energy_particle_interaction_radius);
    
    typedef itk::RCCenterOfMassDirectedFlow<LabelImageType, SegmenterType> COMFlowType;
    COMFlowType::Pointer vCOMFlow = COMFlowType::New();
    vCOMFlow->SetCoefficient(vInternalEnergyCoeffScaler * vParams.energy_coeff_length); // TODO?
    vCOMFlow->SetParticleInteractionRadius(vParams.energy_particle_interaction_radius);
    
    typedef itk::RCShapePrior1DEnergy<LabelImageType, ShapeReferenceImageType> ShapePrior1DType;
    ShapePrior1DType::Pointer vShapePrior1DEnergy = ShapePrior1DType::New();
    vShapePrior1DEnergy->SetBackgroundValue(0);
    vShapePrior1DEnergy->SetCoefficient(vInternalEnergyCoeffScaler * vParams.energy_coeff_shape); 
    vShapePrior1DEnergy->SetTemplateImage(vShapeTemplateReader->GetOutput());

    typedef itk::RCShapePrior23DEnergy<LabelImageType, ShapeReferenceImageType> ShapePrior23DType;
    ShapePrior23DType::Pointer vShapePrior23DEnergy = ShapePrior23DType::New();
    vShapePrior23DEnergy->SetBackgroundValue(0);
    vShapePrior23DEnergy->SetCoefficient(vInternalEnergyCoeffScaler * vParams.energy_coeff_shape_23D); 
    vShapePrior23DEnergy->SetTemplateImage(vShapeTemplateReader->GetOutput());
    
#ifdef USE_GPU
    typedef itk::RCGPUPiecewiseConstantSquareDistCurvatureRegEnergy<LabelImageType, InternalImageType>
    GPUEnergyType;
    GPUEnergyType::Pointer vGPUEnergy = GPUEnergyType::New();
    vGPUEnergy->SetBalloonCoeff(vParams.energy_coeff_balloon);
    vGPUEnergy->SetConstantOutwardFlowCoeff(vParams.energy_coeff_outward_flow);
    vGPUEnergy->SetCurvatureMaskRadius(vParams.energy_curvature_mask_radius);
    vGPUEnergy->SetDataCoeff(vParams.energy_coeff_data);
    vGPUEnergy->SetGPUWorkGroupSize(vParams.gpu_workgroupo_size);
    vGPUEnergy->SetLengthCoeff(vParams.energy_coeff_length);
    vGPUEnergy->SetPiecewiseSmoothMaskRadius(vParams.energy_local_window_radius);
    vGPUEnergy->SetRegionMergingThreshold(vParams.energy_region_merge_ths);
#endif


    /* ENERGIES:
     *      *  Register at segmentation/sampling filters
     */
    SegmenterType::Pointer vSegmentationFilter = SegmenterType::New();
    
    if(vParams.energy_ext_name.compare("pc") == 0) {
        vSegmentationFilter->RegisterExternalEnergy(vPCSqDistEnergy.GetPointer());
    } else if (vParams.energy_ext_name.compare("ps") == 0) {
        vSegmentationFilter->RegisterExternalEnergy(vPSSqDistEnergy.GetPointer());
    } else if (vParams.energy_ext_name.compare("pcDec") == 0) {
        vSegmentationFilter->RegisterExternalEnergy(vPCDecSquaredDistEnergy.GetPointer());
    } else if (vParams.energy_ext_name.compare("pcGauss") == 0) {
        vSegmentationFilter->RegisterExternalEnergy(vPCGaussEnergy.GetPointer());
    } else if (vParams.energy_ext_name.compare("psGauss") == 0) {
        vSegmentationFilter->RegisterExternalEnergy(vPSGaussEnergy.GetPointer());
    } else if (vParams.energy_ext_name.compare("pcDecGauss") == 0) {
        vSegmentationFilter->RegisterExternalEnergy(vPCDecGaussEnergy.GetPointer());
    } else if (vParams.energy_ext_name.compare("pcPoisson") == 0) {
        vSegmentationFilter->RegisterExternalEnergy(vPCPoissonEnergy.GetPointer());
    } else if (vParams.energy_ext_name.compare("psPoisson") == 0) {
        vSegmentationFilter->RegisterExternalEnergy(vPSPoissonEnergy.GetPointer());
    } else if (vParams.energy_ext_name.compare("pcDecPoisson") == 0) {
        vSegmentationFilter->RegisterExternalEnergy(vPCDecPoissonEnergy.GetPointer());
    } else if (vParams.energy_ext_name.compare("psLi") == 0) {
        vSegmentationFilter->RegisterExternalEnergy(vPSLiEnergy.GetPointer());
    } else {
        std::cerr << "Unknown external energy name. Bye bye." << std::endl;
        return 1;
    }
    
    if(vParams.energy_coeff_balloon != 0) {
        vSegmentationFilter->RegisterExternalEnergy(vBalloonFlow.GetPointer());
    }
    
    if(vParams.energy_int_name.compare("manhattenLength") == 0) {
        vSegmentationFilter->RegisterInternalEnergy(vDiscrLengthEnergy.GetPointer());
    } else if(vParams.energy_int_name.compare("curv") == 0) {
        vSegmentationFilter->RegisterInternalEnergy(vKKLengthEnergy.GetPointer());
    } else if(vParams.energy_int_name.compare("springPotential") == 0) {
        vSegmentationFilter->RegisterParticleInteractionEnergy(vSpringPotEnergy.GetPointer());  
    } else if(vParams.energy_int_name.compare("expPotential") == 0) {
        vSegmentationFilter->RegisterParticleInteractionEnergy(vExponentialPotEnergy.GetPointer());  
    } else if(vParams.energy_int_name.compare("centerOfMassDirectedFlow") == 0) {
        vSegmentationFilter->RegisterParticleInteractionEnergy(vCOMFlow.GetPointer());
    } else {
        std::cerr << "Unknown internal energy name. Bye bye." << std::endl;
        return 1;
    }
    if(vParams.energy_coeff_outward_flow != 0) {
        vSegmentationFilter->RegisterInternalEnergy(vConstOutwardFlow.GetPointer());
    }    
    if(vParams.energy_coeff_shape > 0) {
        vSegmentationFilter->RegisterInternalEnergy(vShapePrior1DEnergy.GetPointer());
    }
    if(vParams.energy_coeff_shape_23D) {
        vSegmentationFilter->RegisterInternalEnergy(vShapePrior23DEnergy.GetPointer());    
    }
#ifdef USE_GPU
    // check if the user wants to use the gpu
    if(vParams.seg_use_gpu){ 
        // check if correct energy has been chosen 
        if(vParams.energy_ext_name.compare("ps") == 0 &&
                vParams.energy_int_name.compare("curv") == 0){
            
            vSegmentationFilter->RegisterGPUEnergy(vGPUEnergy.GetPointer());
        } else {
            std::cerr << "Energy not implemented for GPU. " << std::endl;
            return 1;
        }
    }
#endif
    
    SegmenterType::Centroid_Coeff_Type vCentroidCoeff;
    SegmenterType::ArrayType vLocalWindowRadius, vSobolevKernelSigma;
    for (unsigned int vD = 0; vD < DIMENSION; vD++) {
        vLocalWindowRadius[vD] = vParams.energy_local_window_radius  * vSpacingCL[0] * vSpacingCL[0] / vSpacingCL[vD];
        vSobolevKernelSigma[vD] = vParams.energy_sobolev_gradient_kernel_sigma * vSpacingCL[0] * vSpacingCL[0] / vSpacingCL[vD];
        vCentroidCoeff[vD] = vParams.energy_shape_centroid_coeff[vD];
    }

    /*
     * ALGORITHMS:
     *      * Configure
     */
    vSegmentationFilter->SetInput(vDataImagePointer);
    vSegmentationFilter->SetInitImageInput(vInitImagePointer);
    vSegmentationFilter->SetMaxNbIterations(vParams.seg_num_iter);
    vSegmentationFilter->SetAllowFusion(vParams.seg_allow_fusion);
    vSegmentationFilter->SetAllowFission(vParams.seg_allow_fission);
    vSegmentationFilter->SetAllowHandles(vParams.seg_allow_handles);
    vSegmentationFilter->SetUseFastEvolution(vParams.seg_use_fast_evolution);
    vSegmentationFilter->SetUseSobolevGradients(vParams.seg_use_sobolev_gradients);
    vSegmentationFilter->SetSobolevKernelSigma(vParams.energy_sobolev_gradient_kernel_sigma);
    vSegmentationFilter->SetMinimalCellSize(vCellSize);
    vSegmentationFilter->SetDebug(vParams.seg_debug);
    vSegmentationFilter->SetVerbose(vParams.verbose);
    if(vParams.preproc_do_masking) {
        vSegmentationFilter->SetForbiddenRegionInput(vForbiddenImagePointer);
    }


    /*
     * OUTPUT
     *  Rescale and convert such that we get a nice output image
     */
    typedef itk::ShiftScaleImageFilter<LabelAbsImageType, OutputImageType> ScaleForOutputFilterType;

    ScaleForOutputFilterType::Pointer vOutputScaleFilter = ScaleForOutputFilterType::New();
    vOutputScaleFilter->SetScale(1);
    vOutputScaleFilter->SetShift(0);

    vOutputScaleFilter->SetInput(vSegmentationFilter->GetOutput());

    /*
     * OUTPUT:
     *      * Write the image to the output file.
     */
    typedef itk::ImageFileWriter<OutputImageType> FileWriterType;
    FileWriterType::Pointer vFileWriter = FileWriterType::New();
    
    vFileWriter->SetFileName(vOutputFileName);
    vFileWriter->SetInput(vOutputScaleFilter->GetOutput());

    try {
        std::cout << "Labled image is written to: " << vOutputFileName << std::endl;
        vFileWriter->Update();
    } catch (itk::ExceptionObject & e) {
        std::cerr << "Exception caught after starting pipeline in main():" << std::endl;
        std::cerr << e << std::endl;
        std::cerr << "bye bye." << std::endl;
        return 1;
    }

    /*
     * OUTPUT
     *      * Write the parameter file
     */
    vParams.save(vOutputFileName + ".json");
    
    /* 
     * OUTPUT
     *      * End time measurement
     */
    vTimer.Stop();
    std::cout << "Finished! - Time used: " << vTimer.GetTotal() << "s" << std::endl;

    //because --init_mode blob_det doesn't work with --spacing
    //vSpacingCL replaced with vSpacingBadHack
    typedef itk::ChangeInformationImageFilter<LabelAbsImageType> ChangeInformationImageFilterTypeSeg;
    ChangeInformationImageFilterTypeSeg::Pointer vChangeDataImgSpacingFilterSeg = ChangeInformationImageFilterTypeSeg::New();
    vChangeDataImgSpacingFilterSeg->SetInput(vSegmentationFilter->GetOutput());   
    vChangeDataImgSpacingFilterSeg->SetChangeSpacing(true);
    vChangeDataImgSpacingFilterSeg->SetOutputSpacing(vSpacingCL);
    std::cout << "Original image spacing of vSegmentationFilter->GetOutput(): "
              << vSegmentationFilter->GetOutput()->GetSpacing()                     
              << " was overridden to " << vSpacingCL << std::endl;
    vChangeDataImgSpacingFilterSeg->Update();


    typedef itk::ChangeInformationImageFilter<InternalImageType> ChangeInformationImageFilterType23;
    ChangeInformationImageFilterType23::Pointer vChangeDataImgSpacingFilter2 = ChangeInformationImageFilterType23::New();
    InternalImageType::Pointer image3Dlysosomes = extractchannel(image5D, vParams.channels[1]);
    vChangeDataImgSpacingFilter2->SetInput(image3Dlysosomes);   
    vChangeDataImgSpacingFilter2->SetChangeSpacing(true);
    vChangeDataImgSpacingFilter2->SetOutputSpacing(vSpacingCL);
    std::cout << "Original image spacing of image3Dlysosomes: "
              << image3Dlysosomes->GetSpacing()                     
              << " was overridden to " << vSpacingCL << std::endl;
    vChangeDataImgSpacingFilter2->Update();

    typedef itk::LabelImageToShapeLabelMapFilter< LabelAbsImageType, itk::LabelMap
        < itk::StatisticsLabelObject< InternalImageType::PixelType, InternalImageType::ImageDimension > > > 
        LabelImageToShapeLabelMapFilterType;
    typedef itk::StatisticsLabelMapFilter< LabelImageToShapeLabelMapFilterType::OutputImageType, InternalImageType >   
        StatisticsLabelMapFilterType;

    LabelImageToShapeLabelMapFilterType::Pointer labelImageToShapeLabelMapFilter = LabelImageToShapeLabelMapFilterType::New();
    labelImageToShapeLabelMapFilter->SetInput(vChangeDataImgSpacingFilterSeg->GetOutput());
    labelImageToShapeLabelMapFilter->Update();
    StatisticsLabelMapFilterType::Pointer statisticsLabelMapFilter = StatisticsLabelMapFilterType::New();
    statisticsLabelMapFilter->SetInput1(labelImageToShapeLabelMapFilter->GetOutput());
    statisticsLabelMapFilter->SetInput2(vChangeDataImgSpacingFilter2->GetOutput());   //image3Dlysosomes
    statisticsLabelMapFilter->InPlaceOn();
    statisticsLabelMapFilter->Update();

    std::vector<objectStruct> vObjects;
    for(unsigned int i = 0; i < statisticsLabelMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
        objectStruct objectValues;
        StatisticsLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = 
            statisticsLabelMapFilter->GetOutput()->GetNthLabelObject(i);
        objectValues.label = static_cast<int>(labelObject->GetLabel());

        StatisticsLabelMapFilterType::OutputImageType::LabelObjectType::CentroidType centroid = 
            labelObject->GetCentroid ();     
        objectValues.x = centroid[0];
        objectValues.y = centroid[1];
        objectValues.z = centroid[2];
        objectValues.physicalSize = labelObject->GetPhysicalSize();
        objectValues.numberOfPixels = labelObject->GetNumberOfPixels();
        double maxDiameter = 0.0;  
        StatisticsLabelMapFilterType::OutputImageType::LabelObjectType::VectorType ellipsoidVector = 
            labelObject->GetEquivalentEllipsoidDiameter();
        for ( unsigned int vd = 0; vd < ellipsoidVector.GetVectorDimension(); ++vd ) {
            if (maxDiameter < ellipsoidVector[vd]) {
                maxDiameter = ellipsoidVector[vd];
            }
        }
        objectValues.maxDiameter = maxDiameter;
        objectValues.roundness = labelObject->GetRoundness();
        objectValues.mean_lysosome = labelObject->GetMean();
	vObjects.push_back(objectValues);
    }

    //Get mean in Macrophage channel
    ChangeInformationImageFilterType23::Pointer vChangeDataImgSpacingFilter3 = ChangeInformationImageFilterType23::New();
    InternalImageType::Pointer image3Dmacrophage = extractchannel(image5D, vParams.channels[2]);
    vChangeDataImgSpacingFilter3->SetInput(image3Dmacrophage);   
    vChangeDataImgSpacingFilter3->SetChangeSpacing(true);
    vChangeDataImgSpacingFilter3->SetOutputSpacing(vSpacingCL);
    std::cout << "Original image spacing of image3Dmacrophage: "
              << image3Dmacrophage->GetSpacing()                     
              << " was overridden to " << vSpacingCL << std::endl;
    vChangeDataImgSpacingFilter3->Update();

    statisticsLabelMapFilter->SetInput2(vChangeDataImgSpacingFilter3->GetOutput()); 
    statisticsLabelMapFilter->InPlaceOn();
    statisticsLabelMapFilter->Update();

    for(unsigned int i = 0; i < statisticsLabelMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
        StatisticsLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = 
            statisticsLabelMapFilter->GetOutput()->GetNthLabelObject(i);
        vObjects[i].mean_macrophage = labelObject->GetMean();
    }

    boost::filesystem::path p (argv[1]);
    std::string filenameStem = p.stem().string();
    int found;
    found = filenameStem.find_first_of("-");
    std::string coverslipNr = filenameStem.substr(0,found);
    std::string imageNr = filenameStem.substr(found+1);

    //Print values to file
    for(unsigned int i = 0; i < vObjects.size(); ++i) {
        std::ofstream fileout;
        fileout.open(vOutputResultsName.c_str(), std::ofstream::app); 
        fileout << coverslipNr << "\t" << imageNr << "\t" << vObjects[i].label << "\t" 
                << vObjects[i].x << "\t" << vObjects[i].y << "\t" << vObjects[i].z << "\t" 
                << vObjects[i].physicalSize << "\t" << vObjects[i].numberOfPixels << "\t" 
                << vObjects[i].maxDiameter << "\t" << vObjects[i].roundness << "\t"
                << vObjects[i].mean_lysosome << "\t" << vObjects[i].mean_macrophage << "\n";
        fileout.close();
    }


 

    return 0;
}
