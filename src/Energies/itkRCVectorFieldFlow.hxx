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

#include "itkRCVectorFieldFlow.h"

namespace itk
{
    template<typename TLabelImage, typename TDataImage, typename TVectorImage, typename TEnergyDifference>
    RCVectorFieldFlow<TLabelImage, TDataImage, TVectorImage, TEnergyDifference>
    ::RCVectorFieldFlow() {

        m_VecImage = NULL;
    }
    
    template<typename TLabelImage, typename TDataImage, typename TVectorImage, typename TEnergyDifference>
    void RCVectorFieldFlow<TLabelImage, TDataImage, TVectorImage, TEnergyDifference>
    ::PrintSelf(std::ostream & os, Indent indent) const
    {
        Superclass::PrintSelf(os, indent);
        os << indent << "Sum statistic size: " << m_Sums.size() << std::endl;
        os << indent << "Count statistics size: " << m_Count.size() << std::endl;
    }
    
    template<typename TLabelImage, typename TDataImage, typename TVectorImage, typename TEnergyDifference>
    inline 
    typename RCVectorFieldFlow<TLabelImage, TDataImage, TVectorImage, TEnergyDifference>::
    ExternalEnergyReturnType
    RCVectorFieldFlow<TLabelImage, TDataImage, TVectorImage, TEnergyDifference>
    ::EvaluateEnergyDifference_(IndexType aInd, 
            LabelPixelType aLabelFrom, 
            LabelPixelType aLabelTo, 
            DataPixelType aImgValue) {

        unsigned int vNTo = m_Count[aLabelTo];
        unsigned int vNFrom = m_Count[aLabelFrom];
        EnergyDifferenceType vNewToMean = (m_Sums[aLabelTo] + aImgValue) / (vNTo + 1);
        EnergyDifferenceType vMeansFrom = m_Sums[aLabelFrom] / vNFrom;
        EnergyDifferenceType vNewVarTo = this->CalculateVariance(
                m_Sums_2[aLabelTo] + aImgValue * aImgValue, vNewToMean, vNTo+1);
        EnergyDifferenceType vVarFrom = this->CalculateVariance(
                m_Sums_2[aLabelFrom], vMeansFrom, vNFrom);
        

        
        EnergyDifferenceType vEnergy = 0;
        VectorType vFlowVector = m_VecImage->GetPixel(aInd);

        VectorType vNormLabelFrom, vNormLabelTo;
        if (aLabelFrom != 0) { // for the region that is shrinking...
            ApproximateSurfaceNormalFromLabelImg(aInd, aLabelFrom, &vNormLabelFrom);
            vEnergy += vNormLabelFrom * vFlowVector;
        }

        if (aLabelTo != 0) { /// for the region that is growing
            ApproximateSurfaceNormalFromLabelImg(aInd, aLabelTo, &vNormLabelTo);
            vEnergy -= vNormLabelTo * vFlowVector;
        }

        ExternalEnergyReturnType vRetType;

        vRetType.second = this->CalculateKullbackLeiblerDistance(
                vNewToMean,vMeansFrom,vNewVarTo,vVarFrom,vNTo,vNFrom) < 
                Superclass::m_RegionMergingThreshold;
        vRetType.first = vEnergy;
        
        return vRetType;
    }
    
    template<typename TLabelImage, typename TDataImage, typename TVectorImage, typename TEnergyDifference>
    inline void 
    RCVectorFieldFlow<TLabelImage, TDataImage, TVectorImage, TEnergyDifference>
    ::AddPoint_(IndexType aInd, LabelAbsPixelType aLabel, DataPixelType aVal) {

    }
    
    template<typename TLabelImage, typename TDataImage, typename TVectorImage, typename TEnergyDifference>
    inline void 
    RCVectorFieldFlow<TLabelImage, TDataImage, TVectorImage, TEnergyDifference>
    ::RemovePoint_(IndexType aInd /* = 0 */, LabelAbsPixelType aLabel, DataPixelType aVal) {

    }
    
    template<typename TLabelImage, typename TDataImage, typename TVectorImage, typename TEnergyDifference>
    void 
    RCVectorFieldFlow<TLabelImage, TDataImage, TVectorImage, TEnergyDifference>
    ::KillRegion_(LabelAbsPixelType aLabel) {

    }
    
    template<typename TLabelImage, typename TDataImage, typename TVectorImage, typename TEnergyDifference>
    typename RCVectorFieldFlow<TLabelImage, TDataImage, TVectorImage, TEnergyDifference>::EnergyDifferenceType 
    RCVectorFieldFlow<TLabelImage, TDataImage, TVectorImage, TEnergyDifference>
    ::ApproximateSurfaceNormalFromLabelImg(IndexType aIndex, LabelAbsPixelType aL, VectorType* aNormal) {
        
        /// Initilize an neighborhoodoperator for fast access to the label image
        /// in the vicinity of the particle of interest.
        typedef ConstNeighborhoodIterator<LabelImageType> LabelImageNeighborhoodIteratorType;
        typedef typename LabelImageNeighborhoodIteratorType::RadiusType LabelImageNeighborhoodIteratorRadiusType;
        LabelImageNeighborhoodIteratorRadiusType vLabelImageIteratorRadius;
        vLabelImageIteratorRadius.Fill(1);
        LabelImageNeighborhoodIteratorType vLabelImageIterator(
                vLabelImageIteratorRadius,
                this->m_LabelImage,
                this->m_LabelImage->GetBufferedRegion());
        typename LabelImageNeighborhoodIteratorType::NeighborhoodType vNeighhood;

        vLabelImageIterator.SetLocation(aIndex);
        vNeighhood = vLabelImageIterator.GetNeighborhood();

        int vNSize = vNeighhood.Size();

        aNormal->Fill(0);

        for (int vI = 0; vI < vNSize; vI++) {
            OffsetType vOffset = vNeighhood.GetOffset(vI);
            if (static_cast<unsigned int> (abs(vNeighhood.GetElement(vI))) == aL) {
                for (unsigned int vD = 0; vD < LabelImageType::ImageDimension; vD++) {
                    (*aNormal)[vD] -= vOffset[vD];
                }
            }
        }

        EnergyDifferenceType vNorm = aNormal->GetNorm();
        if (vNorm > 0.001) {
            aNormal->Normalize();
        }
        
        
        return vNorm;
    }
        

} // end namespace itk
