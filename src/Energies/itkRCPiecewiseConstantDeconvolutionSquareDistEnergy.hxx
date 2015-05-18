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

#include "itkRCPiecewiseConstantDeconvolutionSquareDistEnergy.h"

namespace itk
{
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    RCPiecewiseConstantDeconvolutionSquareDistEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::RCPiecewiseConstantDeconvolutionSquareDistEnergy() {

    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void RCPiecewiseConstantDeconvolutionSquareDistEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::PrintSelf(std::ostream & os, Indent indent) const
    {
        Superclass::PrintSelf(os, indent);
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>    
    typename RCPiecewiseConstantDeconvolutionSquareDistEnergy<TLabelImage, TDataImage, TEnergyDifference>::
    ExternalEnergyReturnType
    RCPiecewiseConstantDeconvolutionSquareDistEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::EvaluateEnergyDifference_(IndexType aIndex, 
            LabelPixelType aLabelFrom, 
            LabelPixelType aLabelTo, 
            DataPixelType aImgValue) {
        
        
        /*
         * Do not really subtract / add the scaled PSF to the 'model image'.
         * Example for the case where the energyDiff is calculated for a possible
         * change to the background (removal):
         * Calc (J' - I)^2  -  (J - I)^2
         * = (J - (c-BG) * PSF - I)^2 - (J - I)^2
         */

        /// Define the region of influence:

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
        EnergyDifferenceType vEnergyDiff = 0;

        InternalPixelType vIntensity_FromLabel = this->m_Intensities[aLabelFrom];
        InternalPixelType vIntensity_ToLabel = this->m_Intensities[aLabelTo];


        for (vPSFIt.GoToBegin(), vModelIt.GoToBegin(), vDataIt.GoToBegin();
                !vPSFIt.IsAtEnd();
                ++vPSFIt, ++vModelIt, ++vDataIt) {
            InternalPixelType vEOld = (vModelIt.Get() - vDataIt.Get());
            vEOld = vEOld * vEOld;
            //            vEOld = fabs(vEOld);
            InternalPixelType vENew = (vModelIt.Get() - (vIntensity_FromLabel - vIntensity_ToLabel) *
                    vPSFIt.Get()) - vDataIt.Get();
            vENew = vENew * vENew;

            vEnergyDiff += vENew - vEOld;
        }

        ExternalEnergyReturnType vRetType;
        vRetType.first = this->m_Coefficient * vEnergyDiff;
        
        unsigned int vNTo = this->m_Count[aLabelTo];
        unsigned int vNFrom = this->m_Count[aLabelFrom];
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
    


} // end namespace itk
