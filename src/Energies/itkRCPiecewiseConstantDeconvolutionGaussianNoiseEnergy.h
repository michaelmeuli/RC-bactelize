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
#ifndef ITKRCPIECEWISECONSTANTDECONVOLUTIONGAUSSIANNOISEENERGY_H
#define	ITKRCPIECEWISECONSTANTDECONVOLUTIONGAUSSIANNOISEENERGY_H

#include "itkRCPiecewiseConstantDeconvolutionEnergyBaseClass.h"
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
 * constant deconvolution image model and with i.i.d. Gaussian noise. 
 *
 * The energy to minimize is ...
 * 
 */
template<typename TLabelImage, typename TDataImage, typename TEnergyDifference = float>
class ITK_EXPORT RCPiecewiseConstantDeconvolutionGaussianNoiseEnergy :
  public RCPiecewiseConstantDeconvolutionEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
  {
  public :
    /**
     * @name Standard ITK declarations
     */
    //@{
    typedef RCPiecewiseConstantDeconvolutionGaussianNoiseEnergy<TLabelImage, TDataImage, TEnergyDifference> Self;
    typedef RCPiecewiseConstantDeconvolutionEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<Self const> ConstPointer;

    itkNewMacro(Self);
    itkTypeMacro(RCPiecewiseConstantDeconvolutionGaussianNoiseEnergy, 
            RCPiecewiseConstantDeconvolutionBaseClassEnergy);
    
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
    
    inline ExternalEnergyReturnType EvaluateEnergyDifference_(IndexType aInd, 
            LabelPixelType aLabelBefore, 
            LabelPixelType aLabelAfter,
            DataPixelType aImgValue);
   
    /* Protected types */
    protected:
        typedef typename LabelImageType::Pointer LabelImagePointerType;
        typedef typename DataImageType::Pointer DataImagePointerType;
        typedef typename LabelImageType::ConstPointer LabelImageConstPointerType;
        typedef typename DataImageType::ConstPointer DataImageConstPointerType;
        typedef EnergyDifferenceType InternalPixelType;
        typedef Image<InternalPixelType, LabelImageType::ImageDimension> InternalImageType;
        typedef typename InternalImageType::Pointer InternalImagePointerType;
        typedef typename InternalImageType::ConstPointer InternalImageConstPointerType;
        
        /** Statistic typedefs */
        typedef typename Superclass::SumStatisticsMapType SumStatisticsMapType;
        typedef typename SumStatisticsMapType::iterator SumStatisticsIteratorType;
        typedef typename Superclass::CountStatisticsMapType CountStatisticsMapType;
        typedef typename CountStatisticsMapType::iterator CountStatisticsIteratorType;
    
    protected:
        RCPiecewiseConstantDeconvolutionGaussianNoiseEnergy();
        ~RCPiecewiseConstantDeconvolutionGaussianNoiseEnergy() {};
        void PrintSelf(std::ostream & os, Indent indent) const;
        
        /** Reestimates the intensities for the PC regions. This method calls 
         * overwrites but calls the superclass implementation.
         */
        virtual void RenewDeconvolutionStatistics(
                LabelImageConstPointerType aInitImage,
                DataImageConstPointerType aDataImage);
        
    private:
        RCPiecewiseConstantDeconvolutionGaussianNoiseEnergy(const Self &);  //purposely not implemented
        void operator=(const Self &); //purposely not implemented
                

        
        SumStatisticsMapType m_ResidualSum_2;
}; // end class
} // end namespace


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRCPiecewiseConstantDeconvolutionGaussianNoiseEnergy.hxx"
#endif
 
#endif	/* ITKRCPIECEWISECONSTANTDECONVOLUTIONGAUSSIANNOISEENERGY_H */

