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

#include <limits>

#include "itkRCPiecewiseConstantPoissonNoiseEnergy.h"

namespace itk
{
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    RCPiecewiseConstantPoissonNoiseEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::RCPiecewiseConstantPoissonNoiseEnergy() {

    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void RCPiecewiseConstantPoissonNoiseEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::PrintSelf(std::ostream & os, Indent indent) const
    {
        Superclass::PrintSelf(os, indent);
        os << indent << "Sum statistic size: " << this->m_Sums.size() << std::endl;
        os << indent << "Count statistics size: " << this->m_Count.size() << std::endl;
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    inline 
    typename RCPiecewiseConstantPoissonNoiseEnergy<TLabelImage, TDataImage, TEnergyDifference>::
    ExternalEnergyReturnType
    RCPiecewiseConstantPoissonNoiseEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::EvaluateEnergyDifference_(IndexType aInd, 
            LabelPixelType aLabelFrom, 
            LabelPixelType aLabelTo, 
            DataPixelType aImgValue) {
        
        EnergyDifferenceType vNTo = this->m_Count[aLabelTo];
        EnergyDifferenceType vNFrom = this->m_Count[aLabelFrom];
        EnergyDifferenceType vNewToMean = (this->m_Sums[aLabelTo] + aImgValue) / (vNTo + 1);
        EnergyDifferenceType vNewFromMean = (this->m_Sums[aLabelFrom] - aImgValue) / (vNFrom - 1);
        EnergyDifferenceType vOldToMean = (this->m_Sums[aLabelTo]) / (vNTo);
        EnergyDifferenceType vOldFromMean = this->m_Sums[aLabelFrom] / vNFrom;
        
        EnergyDifferenceType vNewToVar = this->CalculateVariance(
                this->m_Sums_2[aLabelTo] + aImgValue * aImgValue, vNewToMean, vNTo+1);
        EnergyDifferenceType vNewFromVar = this->CalculateVariance(
                this->m_Sums_2[aLabelFrom] - aImgValue * aImgValue, vNewFromMean, vNFrom-1);
        EnergyDifferenceType vOldFromVar = this->CalculateVariance(
                this->m_Sums_2[aLabelFrom], vOldFromMean, vNFrom);
        EnergyDifferenceType vOldToVar = this->CalculateVariance(
                this->m_Sums_2[aLabelTo], vOldToMean, vNTo);
        

        if (this->m_Sums_2[aLabelFrom] < 0 || this->m_Sums_2[aLabelTo] < 0) {
            assert(!"Region has negative variance.");
        }
        // it might happen that the variance is exactly 0 (if all pixel contain 
        // the same value). We therefore set it to a small value.
        if(vNewToVar == 0)  vNewToVar = NumericTraits<EnergyDifferenceType>::epsilon()*10;
        if(vOldFromVar == 0) vOldFromVar = NumericTraits<EnergyDifferenceType>::epsilon()*10;
        if(vOldToVar == 0) vOldToVar = NumericTraits<EnergyDifferenceType>::epsilon()*10;
        if(vNewFromVar == 0) vNewFromVar = NumericTraits<EnergyDifferenceType>::epsilon()*10;

        ExternalEnergyReturnType vRetType;
                
        vRetType.first = log(vNewFromMean)*(this->m_Sums[aLabelFrom]-aImgValue)-(vNFrom-1)*vNewFromMean;
        vRetType.first += -(log(vOldFromMean)*(this->m_Sums[aLabelFrom]) - (vNFrom) * vOldFromMean);
        
        vRetType.first += log(vNewToMean)*(this->m_Sums[aLabelTo]+aImgValue)-(vNTo+1)*vNewToMean;
        vRetType.first += -(log(vOldToMean)*(this->m_Sums[aLabelTo]) - (vNTo) * vOldToMean);
        
        vRetType.first *= - this->m_Coefficient;

        
        vRetType.second = this->CalculateKullbackLeiblerDistance(
                vNewToMean,vOldFromMean,vNewToVar,vOldFromVar,vNTo+1,vNFrom) < 
                this->m_RegionMergingThreshold;

        return vRetType;
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    inline void 
    RCPiecewiseConstantPoissonNoiseEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::AddPoint_(IndexType aInd, LabelAbsPixelType aLabel, DataPixelType aVal) {

    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    inline void 
    RCPiecewiseConstantPoissonNoiseEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::RemovePoint_(IndexType aInd /* = 0 */, LabelAbsPixelType aLabel, DataPixelType aVal) {

    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void 
    RCPiecewiseConstantPoissonNoiseEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::KillRegion_(LabelAbsPixelType aLabel) {

    }
        

} // end namespace itk
