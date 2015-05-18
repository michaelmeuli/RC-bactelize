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
#ifndef RCGPUPiecewiseConstantSquareDistCurvatureRegEnergy_H
#define	RCGPUPiecewiseConstantSquareDistCurvatureRegEnergy_H

#include "itkRCGPUEnergyBaseClass.h"
#include "itkConstantPadImageFilter.h"

#ifdef USE_GPU
#include "itkOpenCLSupporter.h"
#endif

#include <map>
namespace itk
{

/**
 * @brief Computes energy differences for a pixel change in for a piece-wise 
 * constant image model and with i.i.d. Gaussian noise. 
 *
 * The energy to minimize is \f$E = \sum_i^M(\mu_i - I(x))^2\f$ with M being the 
 * number of regions, I the image and \f$\mu_i\f$ the mean of region i.
 * 
 * In case of 2 regions this image model is also called Chan-Vese model. 
 * 
 * As for the length term, a curvature-regularizing flow is introduced. 
 * For details, see "J. Kybic and J. Kratky, “Discrete curvature calculation for 
 * fast level set segmentation,” in Proc. 16th IEEE Int. Conf. Image Process., 
 * Nov. 2009, pp. 3017–3020.
 * 
 */
template<typename TLabelImage, typename TDataImage, typename TEnergyDifference = float>
class ITK_EXPORT RCGPUPiecewiseConstantSquareDistCurvatureRegEnergy :
  public RCGPUEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
  {
  public :
    /**
     * @name Standard ITK declarations
     */
    //@{
    typedef RCGPUPiecewiseConstantSquareDistCurvatureRegEnergy<TLabelImage, TDataImage, TEnergyDifference> Self;
    typedef RCGPUEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<Self const> ConstPointer;

    itkNewMacro(Self);
    itkTypeMacro(RCGPUPiecewiseConstantSquareDistCurvatureRegEnergy, RCGPUEnergyBaseClass);
    
    typedef TDataImage DataImageType;
    
    typedef typename Superclass::IndexType IndexType;
    typedef typename Superclass::LabelPixelType LabelPixelType;
    typedef typename Superclass::LabelImageType LabelImageType;
    typedef typename Superclass::LabelAbsPixelType LabelAbsPixelType;
    typedef typename Superclass::DataPixelType DataPixelType;
    typedef typename Superclass::EnergyDifferenceType EnergyDifferenceType;
    typedef typename Superclass::GPU_IntType GPU_IntType;
    typedef typename Superclass::GPU_UIntType GPU_UIntType;
    typedef typename Superclass::SizeType SizeType;

    virtual void EvaluateEnergyDifferences(unsigned int **aGPUCoordinates,
            unsigned int* aGPUToLabel,
            float* aGPUOutInternalEnergy,
            float* aGPUOutExternalEnergy,
            unsigned int* aGPUOutMerge,
            unsigned int aPaddedNbCandidates);
    
    virtual void PrepareEnergyCalculation();
    virtual void PrepareEnergyCalculationForIteration();
    
    /** Overrides base class virtual method. Here: Used to calculate online statistics */
    void AddPoint(IndexType aInd, LabelAbsPixelType aLabel, DataPixelType aVal);
        
    /** Overrides base class virtual method. Here: Used to calculate online statistics . */
    void RemovePoint(IndexType aInd /* = 0 */, LabelAbsPixelType aLabel, DataPixelType aVal);

    /** Implement the virtual method to clean-up the region statistics. */
    void KillRegion(LabelAbsPixelType aL);
    
    void CleanUp();
    
    void SetGPUWorkGroupSize(unsigned int aGPUWorkGroupSize) {
        this->m_gWorkGroupSize = aGPUWorkGroupSize;
    };
    
    itkSetMacro(DataCoeff,float);
    itkGetConstMacro(DataCoeff,float);
    itkSetMacro(LengthCoeff,float);
    itkGetConstMacro(LengthCoeff,float);
    itkSetMacro(BalloonCoeff,float);
    itkGetConstMacro(BalloonCoeff,float);
    itkSetMacro(ConstantOutwardFlowCoeff,float);
    itkGetConstMacro(ConstantOutwardFlowCoeff,float);
    itkSetMacro(RegionMergingThreshold,float);
    itkGetConstMacro(RegionMergingThreshold,float);
    itkSetMacro(CurvatureMaskRadius,float);
    itkGetConstMacro(CurvatureMaskRadius,float);
    itkSetMacro(PiecewiseSmoothMaskRadius,float);
    itkGetConstMacro(PiecewiseSmoothMaskRadius,float);
        
    protected:
        RCGPUPiecewiseConstantSquareDistCurvatureRegEnergy();
        ~RCGPUPiecewiseConstantSquareDistCurvatureRegEnergy() {};
        void PrintSelf(std::ostream & os, Indent indent) const;
        
        
    private:
        RCGPUPiecewiseConstantSquareDistCurvatureRegEnergy(const Self &);  //purposely not implemented
        void operator=(const Self &); //purposely not implemented
            
        /** private methods*/
        void PadLabelImageAndCopyToGPUArray();
        void PadDataImageAndCopyToGPUArray();
        void InitOpenCL();
        void InitGPUMembers();
        

            
        void calculateInternalEnergyWithGPU(
        GPU_IntType* labels, unsigned int* tgsize,
                float* mask_radius, unsigned int* offset,
                float volume, float energy_contour_length_coeff,
                float* curvature_flow, unsigned int* testing,
                unsigned int** candidates, unsigned int* to_label,
                unsigned int size_1D, unsigned int local_size);
             
        void calculateExternalEnergyWithGPU(
        GPU_IntType* labels, float* data,
                unsigned int* tgsize,
                float energy_region_coeff, float region_merging_threshold,
                float balloon_force_coeff, float constant_outward_flow,
                float* mask_radius, unsigned int* offset, float* curvature_flow,
                unsigned int* testing, unsigned int** candidate,
                unsigned int* to_label, unsigned int size_1D, unsigned int local_size,
                float* mean_label, float* var_label, int labels_count);
        
//        void calculateExternalEnergyWithGPU3D(
//                int* labels, float* data, int* masks,
//                unsigned int* tgsize,
//                float energy_region_coeff, float region_merging_threshold,
//                float balloon_force_coeff, float constant_outward_flow,
//                unsigned int* mask_radius, unsigned int* offset, float* curvature_flow,
//                unsigned int* testing, unsigned int* candidate_x, unsigned int* candidate_y, unsigned int* candidate_z,
//                unsigned int* to_label, unsigned int size_1D, unsigned int local_size,
//                float* mean_label, float* var_label, int labels_count);
        
//        EnergyDifferenceType CalculateScaledSphereVolume(float aRadius);
        

        /** parameters */
        float m_DataCoeff;
        float m_LengthCoeff;
        float m_BalloonCoeff;
        float m_ConstantOutwardFlowCoeff;
        float m_RegionMergingThreshold;
        float m_CurvatureMaskRadius;
        float m_PiecewiseSmoothMaskRadius;
        
        /** internals*/
        float m_gCurvatureRadius[LabelImageType::ImageDimension];
        float m_gRegionRadius[LabelImageType::ImageDimension];
        unsigned int m_gWorkGroupSize;

        SizeType m_gDataPaddingSize;
        SizeType m_gLabelPaddingSize;

        GPU_IntType* m_gLabelImage;
        float* m_gDataImage;
        
        cl_program m_gProgram[4];
        cl_kernel m_gKernels[4]; // one for internal and one for external
        cl_platform_id m_gPlatform;
        cl_command_queue m_gCmd_queue;
        cl_context m_gContext;
        cl_device_id m_gDevice;

        
        // energy related  statistics 
        typedef std::map<LabelAbsPixelType, EnergyDifferenceType> SumStatisticsMapType;
        typedef typename SumStatisticsMapType::iterator SumStatisticsIteratorType;
        
        /** Sum of intensities within a region - used to compute means */
        SumStatisticsMapType m_Sums;
        
        /** Sum of squares of intensities within regions - used to compute variances*/
        SumStatisticsMapType m_Sums_2;
        
        typedef std::map<LabelAbsPixelType, unsigned int>  CountStatisticsMapType;
        typedef typename CountStatisticsMapType::iterator CountStatisticsIteratorType;
        
        /** The number of pixel within a regin. */
        CountStatisticsMapType m_Count;
        
}; // end class
} // end namespace


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRCGPUPiecewiseConstantSquareDistCurvatureRegEnergy.hxx"
#endif
 
#endif	/* RCGPUPiecewiseConstantSquareDistCurvatureRegEnergy */

