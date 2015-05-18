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

#include "itkRCPiecewiseConstantDeconvolutionGaussianNoiseEnergy.h"

namespace itk
{
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    RCPiecewiseConstantDeconvolutionGaussianNoiseEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::RCPiecewiseConstantDeconvolutionGaussianNoiseEnergy() {

    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void RCPiecewiseConstantDeconvolutionGaussianNoiseEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::PrintSelf(std::ostream & os, Indent indent) const
    {
        Superclass::PrintSelf(os, indent);
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>    
    typename RCPiecewiseConstantDeconvolutionGaussianNoiseEnergy<TLabelImage, TDataImage, TEnergyDifference>::
    ExternalEnergyReturnType
    RCPiecewiseConstantDeconvolutionGaussianNoiseEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::EvaluateEnergyDifference_(IndexType aIndex, 
            LabelPixelType aLabelFrom, 
            LabelPixelType aLabelTo, 
            DataPixelType aImgValue) {
        

        RegionType vRegion;
        SizeType vRegionSize = this->m_PSF->GetLargestPossibleRegion().GetSize();

        OffsetType vOffset;
        vOffset.Fill(vRegionSize[0] / 2); // we assume the region to be a hypercube.
        vRegion.SetSize(vRegionSize);
        vRegion.SetIndex(aIndex - vOffset);

        /// Move the PSF window such that it is centered on aIndex
        RegionType vPSFRegion = this->m_PSF->GetLargestPossibleRegion().GetSize();
        vPSFRegion.SetIndex(vRegion.GetIndex());
        this->m_PSF->SetBufferedRegion(vPSFRegion);

        /// After the cropping at the data-image boundaries, vRegion will be the
        /// region to treat in the data-image space.
        vRegion.Crop(this->m_DeconvolutionModelImage->GetBufferedRegion());

        /// Iterate through the region and subtract the psf from the conv image.
        ImageRegionConstIterator<InternalImageType> vModelIt(this->m_DeconvolutionModelImage, vRegion);
        ImageRegionConstIterator<InternalImageType> vPSFIt(this->m_PSF, vRegion);
        ImageRegionConstIterator<InternalImageType> vDataIt(this->m_DataImage, vRegion);
        ImageRegionConstIterator<LabelImageType> vLabelIt(this->m_LabelImage, vRegion);


        InternalPixelType vIntensity_FromLabel = this->m_Intensities[aLabelFrom];
        InternalPixelType vIntensity_ToLabel = this->m_Intensities[aLabelTo];

        EnergyDifferenceType vResidualAtX = (aImgValue - this->m_DeconvolutionModelImage->GetPixel(aIndex));
        EnergyDifferenceType vOneBySq2Pi = 1.0/sqrt(2.0*M_PI);
        EnergyDifferenceType vNFrom = this->m_Count[aLabelFrom];
        EnergyDifferenceType vNTo = this->m_Count[aLabelTo];

        std::set<LabelAbsPixelType> vLabelsInSupport;

        SumStatisticsMapType vResidualSum2Old;
        SumStatisticsMapType vResidualSum2New;
//        SumStatisticsMapType vN;
        for (vPSFIt.GoToBegin(), vModelIt.GoToBegin(), vDataIt.GoToBegin(), vLabelIt.GoToBegin();
                !vPSFIt.IsAtEnd();
                ++vPSFIt, ++vModelIt, ++vDataIt, ++vLabelIt) {

            LabelAbsPixelType vAbsLabel = static_cast<LabelAbsPixelType>(abs(vLabelIt.Get()));
            if(vLabelsInSupport.find(vAbsLabel)==vLabelsInSupport.end()) {
                vLabelsInSupport.insert(vAbsLabel);
                vResidualSum2Old[vAbsLabel] = 0;
                vResidualSum2New[vAbsLabel] = 0;
            }
            
            vResidualSum2Old[vAbsLabel] += (vModelIt.Get() - vDataIt.Get()) * (vModelIt.Get() - vDataIt.Get());

            EnergyDifferenceType vRNew = (vModelIt.Get() + (vIntensity_ToLabel - vIntensity_FromLabel) *
                    vPSFIt.Get()) - vDataIt.Get();
            vResidualSum2New[vAbsLabel] += vRNew * vRNew;

//            vN[vAbsLabel]++;
        }

        EnergyDifferenceType vEnergyDiff = 0;
        for(typename std::set<LabelAbsPixelType>::iterator vLabelsIt = vLabelsInSupport.begin();
                vLabelsIt != vLabelsInSupport.end(); ++vLabelsIt) {
  
            EnergyDifferenceType vNold = this->m_Count[*vLabelsIt];
            EnergyDifferenceType vNnew = vNold;
            if(*vLabelsIt == aLabelFrom) {
                vNnew--;
            }else if(*vLabelsIt == aLabelTo) {
                vNnew++;
            }
                
            EnergyDifferenceType vSigmaNew = sqrt(
                    (m_ResidualSum_2[*vLabelsIt] - vResidualSum2Old[*vLabelsIt] 
                    + vResidualSum2New[*vLabelsIt]) / (vNnew-1));
            EnergyDifferenceType vSigmaOld = sqrt(
                    (m_ResidualSum_2[*vLabelsIt]) / (vNold-1));
            
            vEnergyDiff += vNnew*log(vOneBySq2Pi/vSigmaNew);
            vEnergyDiff += -vNold*log(vOneBySq2Pi/vSigmaOld);
            vEnergyDiff += vResidualSum2New[*vLabelsIt] / (2 * vSigmaNew*vSigmaNew);
            vEnergyDiff += -vResidualSum2Old[*vLabelsIt] / (2 * vSigmaOld*vSigmaOld);
        }
        vEnergyDiff *= -1.0f;

        
        ExternalEnergyReturnType vRetType;
        vRetType.first = this->m_Coefficient * vEnergyDiff;
        
//        unsigned int vNTo = this->m_Count[aLabelTo];
//        unsigned int vNFrom = this->m_Count[aLabelFrom];
        EnergyDifferenceType vNewToMean = (this->m_Sums[aLabelTo] + aImgValue) / (vNTo + 1);
        EnergyDifferenceType vMeansFrom = this->m_Sums[aLabelFrom] / vNFrom;
        EnergyDifferenceType vNewVarTo = this->CalculateVariance(
                this->m_Sums_2[aLabelTo] + aImgValue * aImgValue, vNewToMean, vNTo+1);
        EnergyDifferenceType vVarFrom = this->CalculateVariance(
                this->m_Sums_2[aLabelFrom], vMeansFrom, vNFrom);
        
        vRetType.second = this->CalculateKullbackLeiblerDistance(
                vNewToMean,vMeansFrom,vNewVarTo,vVarFrom,vNTo,vNFrom) < 
                this->m_RegionMergingThreshold;
        
        return vRetType;
    }
    
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void RCPiecewiseConstantDeconvolutionGaussianNoiseEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::RenewDeconvolutionStatistics(LabelImageConstPointerType aInitImage, 
            DataImageConstPointerType aDataImage) {
        
        Superclass::RenewDeconvolutionStatistics(aInitImage, aDataImage);
        
        // The superclass method has renewed intensities and the model image J.
        // Here we additionally need the variance of (I-J) region wise.
        
        // Reinitalize the residual sum of squares
        m_ResidualSum_2.clear();

        for(CountStatisticsIteratorType vIt = this->m_Count.begin(); 
                vIt != this->m_Count.end();++vIt) {
            m_ResidualSum_2[vIt->first] = 0;
        }
        
        ImageRegionConstIterator<InternalImageType> vModelIt(this->m_DeconvolutionModelImage, 
                this->m_DeconvolutionModelImage->GetBufferedRegion());
        ImageRegionConstIterator<InternalImageType> vDataIt(this->m_DataImage, 
                this->m_DataImage->GetBufferedRegion());
        ImageRegionConstIterator<LabelImageType> vLabelIt(this->m_LabelImage, 
                this->m_LabelImage->GetBufferedRegion());
                
        vModelIt.GoToBegin();
        vDataIt.GoToBegin();
        vLabelIt.GoToBegin();
        for (; !vModelIt.IsAtEnd(); ++vModelIt, ++vDataIt, ++vLabelIt) {
            m_ResidualSum_2[abs(vLabelIt.Get())] += (vDataIt.Get() - vModelIt.Get()) * (vDataIt.Get() - vModelIt.Get());
        }
        
    }

} // end namespace itk
