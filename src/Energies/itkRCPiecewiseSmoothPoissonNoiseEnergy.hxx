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

#include "itkRCPiecewiseSmoothPoissonNoiseEnergy.h"

namespace itk
{
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    RCPiecewiseSmoothPoissonNoiseEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::RCPiecewiseSmoothPoissonNoiseEnergy() {
        m_Radius = 10;
        m_SphereMaskForLocalEnergy = SphereMaskImageSourceType::New();
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void RCPiecewiseSmoothPoissonNoiseEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::PrintSelf(std::ostream & os, Indent indent) const
    {
        Superclass::PrintSelf(os, indent);
        os << indent << "Radius: " << m_Radius << std::endl;
    }
    
        
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void RCPiecewiseSmoothPoissonNoiseEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::PrepareEnergyCalculation() {
        typename SphereMaskImageSourceType::SizeType vMaskSize;
        typename SphereMaskImageSourceType::ArrayType vRadius;
        
        for (unsigned int vD = 0; vD < LabelImageType::ImageDimension; vD++) {
            // the radius is expected to be given in px size of the first 
            // axis. We scale for the all the following dimensions according
            // to the image spacing.
            float vScaledRadius = m_Radius *
            this->m_DataImage->GetSpacing()[0] /
            this->m_DataImage->GetSpacing()[vD];
            
            vMaskSize[vD] = 2 * ceil(vScaledRadius) + 1;
            vRadius[vD] = vScaledRadius;
        }
        
        m_SphereMaskForLocalEnergy->SetRadius(vRadius);
        m_SphereMaskForLocalEnergy->SetSize(vMaskSize);
        m_SphereMaskForLocalEnergy->Update();
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    inline 
    typename RCPiecewiseSmoothPoissonNoiseEnergy<TLabelImage, TDataImage, TEnergyDifference>::
    ExternalEnergyReturnType
    RCPiecewiseSmoothPoissonNoiseEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::EvaluateEnergyDifference_(IndexType aInd, 
            LabelPixelType aLabelFrom, 
            LabelPixelType aLabelTo, 
            DataPixelType aImgValue) {
        
        // vRegion is the size of our temporary window
        typedef typename LabelImageType::RegionType RegionType;
        RegionType vRegion;
        
        // vOffset is basically the difference of the center and the start of the window
        typedef typename LabelImageType::OffsetType OffsetType;
        OffsetType vOffset;

        //read out the size of the mask
        vRegion = m_SphereMaskForLocalEnergy->GetOutput()->GetBufferedRegion();
        for (unsigned int vD = 0; vD < LabelImageType::ImageDimension; vD++) {
            // flooring is on purpose!
            vOffset[vD] = (vRegion.GetSize())[vD] / 2;
        }
        // translate the region
        vRegion.SetIndex(aInd - vOffset);
        // apply translation to the mask
        m_SphereMaskForLocalEnergy->GetOutput()->SetBufferedRegion(vRegion);
        // Crop the region with the data domain
        vRegion.Crop(this->m_DataImage->GetBufferedRegion());
        // we assume the label image to be in memory and to be as large as the
        // data image. The next line is most probably not necessary. //TODO
        vRegion.Crop(this->m_LabelImage->GetBufferedRegion());


        ImageRegionConstIterator<LabelImageType> vLabelIt(this->m_LabelImage, vRegion);
        ImageRegionConstIterator<DataImageType> vDataIt(this->m_DataImage, vRegion);
        ImageRegionConstIterator<MaskImageType> vMaskIt(
                m_SphereMaskForLocalEnergy->GetOutput(),
                vRegion); //m_SphereMaskForLocalEnergy->GetOutput()->GetLargestPossibleRegion());

        EnergyDifferenceType vSumFrom = -aImgValue; // we ignore the value of the center point
        EnergyDifferenceType vSumTo = 0;
        EnergyDifferenceType vSumOfSqFrom = -aImgValue*aImgValue; // ignore the value of the center point.
        EnergyDifferenceType vSumOfSqTo = 0.0;
        EnergyDifferenceType vNFrom = -1.0;
        EnergyDifferenceType vNTo = 0;

        for (vLabelIt.GoToBegin(), vDataIt.GoToBegin(), vMaskIt.GoToBegin();
                !vLabelIt.IsAtEnd();
                ++vLabelIt, ++vDataIt, ++vMaskIt) {

            if (vMaskIt.Get()) {
                if (static_cast<unsigned int> (abs(vLabelIt.Get())) == aLabelFrom) {
                    vSumFrom += vDataIt.Get();
                    vSumOfSqFrom += vDataIt.Get() * vDataIt.Get();
                    vNFrom++;
                } else if (static_cast<unsigned int> (abs(vLabelIt.Get())) == aLabelTo) {
                    vSumTo += vDataIt.Get();
                    vSumOfSqTo += vDataIt.Get() * vDataIt.Get();
                    vNTo++;
                }
            }
        }
        EnergyDifferenceType vMeanTo;
        EnergyDifferenceType vVarTo;
        EnergyDifferenceType vMeanFrom;
        EnergyDifferenceType vVarFrom;
        if (vNTo == 0) { // this should only happen with the BG label
            EnergyDifferenceType vDivisor = this->m_Count[aLabelTo];
            if(vDivisor == 0) vDivisor = 1;
            vMeanTo = this->m_Sums[aLabelTo]/vDivisor;
            vVarTo = this->CalculateVariance(this->m_Sums_2[aLabelTo], vMeanTo,vNTo);
        } else {
            vMeanTo = vSumTo / vNTo;
            //            vVarTo = (vSumOfSqTo - vSumTo * vSumTo / vNTo) / (vNTo - 1);
            vVarTo = this->CalculateVariance(vSumOfSqTo, vMeanTo, vNTo);
        }

        if (vNFrom == 0) {
            EnergyDifferenceType vDivisor = this->m_Count[aLabelFrom];
            if(vDivisor == 0) vDivisor = 1;
            vMeanFrom = this->m_Sums[aLabelFrom]/vDivisor;
            vVarFrom = this->CalculateVariance(this->m_Sums_2[aLabelFrom],vMeanFrom,vNFrom);
        } else {
            vMeanFrom = vSumFrom / vNFrom;
            vVarFrom = this->CalculateVariance(vSumOfSqFrom, vMeanFrom, vNFrom);
        }

        ExternalEnergyReturnType vReturnVal;
        vReturnVal.second = false;
        if (aLabelFrom != 0 && aLabelTo != 0) {
            if (this->CalculateKullbackLeiblerDistance(vMeanFrom, vMeanTo,
                    vVarFrom, vVarTo, vNFrom, vNTo) < this->m_RegionMergingThreshold) {
                vReturnVal.second = true;
            }
        }
        
        EnergyDifferenceType vEnergyDiff = 0; 

        if (vMeanFrom < NumericTraits<EnergyDifferenceType>::epsilon() &&
                vMeanTo < NumericTraits<EnergyDifferenceType>::epsilon()) {
            vEnergyDiff = 0;
        } else if (vMeanFrom <= NumericTraits<EnergyDifferenceType>::epsilon()) {
            vEnergyDiff = NumericTraits<EnergyDifferenceType>::min();
        } else if (vMeanTo <= NumericTraits<EnergyDifferenceType>::epsilon()) {
            vEnergyDiff = NumericTraits<EnergyDifferenceType>::max();
        } else {
            vEnergyDiff = (aImgValue) * log((vMeanFrom) / (vMeanTo))
                    - (vMeanFrom) + (vMeanTo);
        }
        
        vReturnVal.first = this->m_Coefficient * (vEnergyDiff);

        return vReturnVal;
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    inline void 
    RCPiecewiseSmoothPoissonNoiseEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::AddPoint_(IndexType aInd, LabelAbsPixelType aLabel, DataPixelType aVal) {

    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    inline void 
    RCPiecewiseSmoothPoissonNoiseEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::RemovePoint_(IndexType aInd /* = 0 */, LabelAbsPixelType aLabel, DataPixelType aVal) {
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void 
    RCPiecewiseSmoothPoissonNoiseEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::KillRegion_(LabelAbsPixelType aLabel) {
      
    }
        

} // end namespace itk