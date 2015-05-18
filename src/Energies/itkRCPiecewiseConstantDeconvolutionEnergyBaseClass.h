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
#ifndef ITKRCPIECEWISECONSTANTDECONVOLUTIONENERGYBASECLASS_H
#define	ITKRCPIECEWISECONSTANTDECONVOLUTIONENERGYBASECLASS_H

#include "itkRCExternalEnergyBaseClass.h"
#include "itkFFTConvolutionImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkDivideImageFilter.h"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/stats.hpp>

#include <map>


namespace itk
{

/**
 * @brief Computes energy differences for a pixel change in for a piece-wise 
 * constant deconvolution image model.
 *
 * The energy to minimize is ...
 * 
 */
template<typename TLabelImage, typename TDataImage, typename TEnergyDifference = float>
class ITK_EXPORT RCPiecewiseConstantDeconvolutionEnergyBaseClass :
  public RCExternalEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
  {
  public :
    /**
     * @name Standard ITK declarations
     */
    //@{
    typedef RCPiecewiseConstantDeconvolutionEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference> Self;
    typedef RCExternalEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<Self const> ConstPointer;

    //    itkNewMacro(Self);
    itkTypeMacro(RCPiecewiseConstantDeconvolutionEnergyBaseClass, RCExternalEnergyBaseClass);
    
    typedef typename Superclass::IndexType IndexType;
    typedef typename Superclass::LabelPixelType LabelPixelType;
    typedef typename Superclass::LabelAbsPixelType LabelAbsPixelType;
    typedef typename Superclass::DataPixelType DataPixelType;
    typedef typename Superclass::EnergyDifferenceType EnergyDifferenceType;
    typedef typename Superclass::ExternalEnergyReturnType ExternalEnergyReturnType;
    typedef typename Superclass::RegionType RegionType;
    typedef typename Superclass::SizeType SizeType; 
    typedef typename Superclass::OffsetType OffsetType;
    typedef typename Superclass::DataImageType DataImageType;
    typedef typename Superclass::LabelImageType LabelImageType;
    typedef typename Superclass::ReconstructedImageType ReconstructedImageType;
    
    ExternalEnergyReturnType EvaluateEnergyDifference_(IndexType aInd, 
            LabelPixelType aLabelBefore, 
            LabelPixelType aLabelAfter,
            DataPixelType aImgValue) = 0;
   
    /** Add a point to the statistics of region. */
    inline void AddPoint_(IndexType aInd, LabelAbsPixelType aRegion, DataPixelType aVal);
    
    /** Remove a point from the statistics of a region. */
    inline void RemovePoint_(IndexType aInd, LabelAbsPixelType aRegion, DataPixelType aVal);
    
    /** Overwrites the implementation of the base class. */
    inline void SwitchPoint(IndexType aInd, LabelAbsPixelType aLabelFrom, 
            LabelAbsPixelType aLabelTo, DataPixelType aVal);

    /** Implement the virtual method to clean-up the region statistics. */
    void KillRegion_(LabelAbsPixelType aL);
    
    /** Implement the virtual method */
    void PrepareEnergyCalculation();
    
    /** Implement the virtual method */
    void PrepareEnergyCalculationForIteration();
    
    itkSetMacro(DoNormalizePSF, bool);
    itkGetMacro(DoNormalizePSF, bool);
    itkSetMacro(OptimizationMode, bool);
    itkGetMacro(OptimizationMode, bool);
    
    /** Get the PSF image. */
    const DataImageType * GetDataImage() const
    { return m_PSFInput.GetPointer(); }
    
    /** Set the label image. This needs to be called to define a valid 
     *  external energy.
     */
    void SetPSF(const DataImageType *aPtr)
    { m_PSFInput = aPtr; }

    /** override the base class method*/
    EnergyDifferenceType CalculateTotalEnergy();
    
    /** Overrides the base class method */
    void GenerateReconstructedImage(void* aPointerToResultImage);

    
    protected:
        RCPiecewiseConstantDeconvolutionEnergyBaseClass();
        ~RCPiecewiseConstantDeconvolutionEnergyBaseClass() {};
        void PrintSelf(std::ostream & os, Indent indent) const;
        
        inline void AddPoint(IndexType aInd, LabelAbsPixelType aRegion, 
                DataPixelType aVal, bool aDoUpdateModelImage);
        inline void RemovePoint(IndexType aInd, LabelAbsPixelType aRegion, DataPixelType aVal,
                bool aDoUpdateModelImage);
    private:
        RCPiecewiseConstantDeconvolutionEnergyBaseClass(const Self &);  //purposely not implemented
        void operator=(const Self &); //purposely not implemented
           
   protected:
        /** Protected typedefs */
        typedef typename LabelImageType::Pointer LabelImagePointerType;
        typedef typename DataImageType::Pointer DataImagePointerType;
        typedef typename LabelImageType::ConstPointer LabelImageConstPointerType;
        typedef typename DataImageType::ConstPointer DataImageConstPointerType;
        typedef EnergyDifferenceType InternalPixelType;
        typedef Image<InternalPixelType, LabelImageType::ImageDimension> InternalImageType;
        typedef typename InternalImageType::Pointer InternalImagePointerType;
        typedef typename InternalImageType::ConstPointer InternalImageConstPointerType;
        
        /** Statistic typedefs */
        typedef std::map<LabelAbsPixelType, EnergyDifferenceType> SumStatisticsMapType;
        typedef typename SumStatisticsMapType::iterator SumStatisticsIteratorType;
        typedef std::map<LabelAbsPixelType, unsigned int>  CountStatisticsMapType;
        typedef typename CountStatisticsMapType::iterator CountStatisticsIteratorType;
        
        /** Private images */
        InternalImagePointerType m_DeconvolutionModelImage;
        InternalImagePointerType m_PSF;
        DataImageConstPointerType m_PSFInput;
        bool m_DoNormalizePSF;
        
        /** Reestimates the intensities for the PC regions. */
        virtual void RenewDeconvolutionStatistics(
                LabelImageConstPointerType aInitImage,
                DataImageConstPointerType aDataImage);
        
        /** Generate the model image J. */
        void GenerateModelImage(
                void* aPointerToResultImage,
                LabelImageConstPointerType aLabelImage,
                InternalImageConstPointerType aDataImage,
                InternalImagePointerType aPSFImage,
                SumStatisticsMapType* aIntensityValues);

        
        /** Update the model image when a point is added, removed or switched. */
        void UpdateModelImage(IndexType aIndex,
                LabelAbsPixelType aFromLabel,
                LabelAbsPixelType aToLabel);
        
        
        /** Estimated intensity values for PC regions - used to generate J */
        SumStatisticsMapType m_Intensities;
        

        
        protected:
            bool m_SwitchPointMode;
            bool m_OptimizationMode; // or sampling mode

            typedef std::set<LabelAbsPixelType> LabelAbsSetType;
            typedef typename LabelAbsSetType::iterator LabelAbsSetIteratorType;
            LabelAbsSetType m_NewRegions;
}; // end class
} // end namespace


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRCPiecewiseConstantDeconvolutionEnergyBaseClass.hxx"
#endif
 
#endif	/* ITKRCPIECEWISECONSTANTDECONVOLUTIONENERGYBASECLASS_H */

