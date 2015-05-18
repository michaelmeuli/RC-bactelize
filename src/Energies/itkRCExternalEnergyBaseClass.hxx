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
#include "itkRCExternalEnergyBaseClass.h"
#include "vnl/vnl_math.h"

namespace itk
{
    template <class TLabelImage, class TDataImage, class TEnergyDifference>
    RCExternalEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
    ::RCExternalEnergyBaseClass() {//: RCEnergyBaseClass() {
        m_DataImage = NULL;
        m_RegionMergingThreshold = 1;
    }
    
    template <class TLabelImage, class TDataImage, class TEnergyDifference>
    void RCExternalEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
    ::PrintSelf(std::ostream & os, Indent indent) const
    {
        Superclass::PrintSelf(os, indent);
        os << indent << "Region merging threshold: " << m_RegionMergingThreshold << std::endl;
//        os << indent << "Data image ptr: " << m_DataImage << std::endl;
    }
    
    template <class TLabelImage, class TDataImage, class TEnergyDifference>
    typename RCExternalEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
    ::ExternalEnergyReturnType
    RCExternalEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
    ::EvaluateEnergyDifference(IndexType aInd, 
                LabelPixelType aLabelBefore, 
                LabelPixelType aLabelAfter,
                DataPixelType aImgValue){
        
        /// The energy is defined to be -\inf in case the particle is the last
        /// parent:
        if(m_Count[aLabelBefore] == 1) {
            /// This works for symmetric types (?)
            return ExternalEnergyReturnType(
                    -NumericTraits<EnergyDifferenceType>::max(), false);
        }
        if(m_Count[aLabelAfter] == 0) { // dangling particles
            return ExternalEnergyReturnType(
                    NumericTraits<EnergyDifferenceType>::max(), false);
        }

        return EvaluateEnergyDifference_(aInd, aLabelBefore, aLabelAfter, aImgValue);
    }
    
    template <class TLabelImage, class TDataImage, class TEnergyDifference>
    void RCExternalEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
    ::EvaluateEnergyDifferences(const std::vector<IndexType>& aInds, 
            const std::vector<LabelPixelType>& aLabelBefore, 
            const std::vector<LabelPixelType>& aLabelAfter,
            const std::vector<DataPixelType>& aImgValue,
            std::vector<ExternalEnergyReturnType>& aReturnVec) {
        
        aReturnVec.clear();
        aReturnVec.resize(aInds.size());
        
        for(unsigned int vI = 0; vI < aReturnVec.size(); vI++) {
            aReturnVec[vI] =
                    EvaluateEnergyDifference_(aInds[vI], 
                    aLabelBefore[vI], 
                    aLabelAfter[vI], 
                    aImgValue[vI]);
            
        }
    }
    
    template <class TLabelImage, class TDataImage, class TEnergyDifference>
    inline  
    typename RCExternalEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
    ::EnergyDifferenceType 
    RCExternalEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
    ::CalculateKullbackLeiblerDistance( EnergyDifferenceType aMu1, EnergyDifferenceType aMu2, 
            EnergyDifferenceType aVar1, EnergyDifferenceType aVar2,
            EnergyDifferenceType aN1, EnergyDifferenceType aN2) {

        EnergyDifferenceType vMu12 = (aN1 * aMu1 + aN2 * aMu2) /
                (aN1 + aN2);

        EnergyDifferenceType vSumOfSq1 = aVar1 * (aN1 - 1) + aN1 * aMu1 * aMu1;
        EnergyDifferenceType vSumOfSq2 = aVar2 * (aN2 - 1) + aN2 * aMu2 * aMu2;

        EnergyDifferenceType vVar12 = (1.0 / (aN1 + aN2 - 1.0)) *
                (vSumOfSq1 + vSumOfSq2 - (aN1 + aN2) * vMu12 * vMu12);

        EnergyDifferenceType vDKL1 = (aMu1 - vMu12) * (aMu1 - vMu12) / (2.0f * vVar12)
                + 0.5f * (aVar1 / vVar12 - 1.0f - log(aVar1 / vVar12));

        EnergyDifferenceType vDKL2 = (aMu2 - vMu12) * (aMu2 - vMu12) / (2.0f * vVar12)
                + 0.5f * (aVar2 / vVar12 - 1.0f - log(aVar2 / vVar12));

//        if(NumericTraits<EnergyDifferenceType>::quiet_NaN())
        // TODO: handle divisions by zero: is this meaningful?
        if(vnl_math_isnan(vDKL1 + vDKL2)) {
            return NumericTraits<EnergyDifferenceType>::max();
        }
        return vDKL1 + vDKL2;
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void 
    RCExternalEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
    ::AddPoint(IndexType aInd, LabelAbsPixelType aLabel, DataPixelType aVal) {
        CountStatisticsIteratorType vIt = m_Count.find(aLabel);
        if(vIt != m_Count.end()) {
            vIt->second++;
            m_Sums[aLabel] += aVal;
            m_Sums_2[aLabel] += aVal * aVal;
        } else {
            m_Count[aLabel] = 1;
            m_Sums[aLabel] = aVal;
            m_Sums_2[aLabel] = aVal * aVal;
        }
        AddPoint_(aInd, aLabel, aVal);
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    inline void 
    RCExternalEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
    ::RemovePoint(IndexType aInd /* = 0 */, LabelAbsPixelType aLabel, DataPixelType aVal) {
        assert(m_Count.find(aLabel) != m_Count.end());
        if(--m_Count[aLabel] == 0) {
            KillRegion_(aLabel);
        } else {
            m_Sums[aLabel] -= aVal;
            m_Sums_2[aLabel] -= aVal * aVal;
        }
        RemovePoint_(aInd, aLabel, aVal);
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void 
    RCExternalEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
    ::KillRegion(LabelAbsPixelType aLabel) {
        if(aLabel != 0) {
            m_Count.erase(aLabel);
            m_Sums.erase(aLabel);
            m_Sums_2.erase(aLabel);
        } else {
            m_Count[0] = 0;
            m_Sums[0] = 0;
            m_Sums_2[0] = 0;
        }
        KillRegion_(aLabel);
    }
    
    

} // end namespace itk
