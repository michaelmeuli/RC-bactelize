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

#include "itkRCDiscreteContourLengthApproxEnergy.h"
#include "itkRCEnergyBaseClass.h"


namespace itk
{
    template<typename TLabelImage, typename TEnergyDifference>
    RCDiscreteContourLengthApproxEnergyEnergy<TLabelImage, TEnergyDifference>
    ::RCDiscreteContourLengthApproxEnergyEnergy() {
        
        ConstNeighborhoodIterator<LabelImageType> vImgIt;
        typename LabelImageType::SizeType vSize;

        if(LabelImageType::ImageDimension == 2) {

            m_lengthWeightMask.resize(16);
            m_Offsets.resize(16);

            // m_Offsets[0] = {{0,1}};
            // m_Offsets[1] = {{0,-1}};
            // m_Offsets[2] = {{1,0}};
            // m_Offsets[3] = {{-1,0}};
            // m_Offsets[4] = {{1,1}};
            // m_Offsets[5] = {{1,-1}};
            // m_Offsets[6] = {{-1,1}};
            // m_Offsets[7] = {{-1,-1}};
            // m_Offsets[8] = {{1,2}};
            // m_Offsets[9] = {{1,-2}};
            // m_Offsets[10] = {{-1,2}};
            // m_Offsets[11] = {{-1,-2}};
            // m_Offsets[12] = {{2,1}};
            // m_Offsets[13] = {{2,-1}};
            // m_Offsets[14] = {{-2,1}};
            // m_Offsets[15] = {{-2,-1}};
            
            m_Offsets[0][0] = 0;
            m_Offsets[0][1] = 1;
            m_Offsets[1][0] = 0;
            m_Offsets[1][1] = -1;
            m_Offsets[2][0] = 1;
            m_Offsets[2][1] = 0;
            m_Offsets[3][0] = -1;
            m_Offsets[3][1] = 0;
            m_Offsets[4][0] = 1;
            m_Offsets[4][1] = 1;
            m_Offsets[5][0] = 1;
            m_Offsets[5][1] = -1;
            m_Offsets[6][0] = -1;
            m_Offsets[6][1] = 1;
            m_Offsets[7][0] = -1;
            m_Offsets[7][1] = -1;
            m_Offsets[8][0] = 1;
            m_Offsets[8][1] = 2;
            m_Offsets[9][0] = 1;
            m_Offsets[9][1] = -2;
            m_Offsets[10][0] = -1;
            m_Offsets[10][1] = 2;
            m_Offsets[11][0] = -1;
            m_Offsets[11][1] = -2;
            m_Offsets[12][0] = 2;
            m_Offsets[12][1] = 1;
            m_Offsets[13][0] = 2;
            m_Offsets[13][1] = -1;
            m_Offsets[14][0] = -2;
            m_Offsets[14][1] = 1;
            m_Offsets[15][0] = -2;
            m_Offsets[15][1] = -1;

            vSize.Fill(2);
            vImgIt.SetRadius(vSize);
            m_Radius = 2;
        } else {
            
            m_lengthWeightMask.resize(static_cast<unsigned int>(pow(LabelImageType::ImageDimension,3)));
            m_Offsets.resize(static_cast<unsigned int>(pow(LabelImageType::ImageDimension,3)));
            vSize.Fill(1);            
            vImgIt.SetRadius(vSize); 

            for(int vI = 0;vI <  vImgIt.Size();vI++) {
                m_Offsets[vI] = vImgIt.GetOffset(vI);
            }
            m_Radius = 1;
        }
        
        
        m_lengthWeightMask.resize(m_Offsets.size());
        for (unsigned int vN = 0; vN < m_Offsets.size(); vN++) {
            OffsetType vOff = m_Offsets[vN];
            float vLength = 0.0f;
            for(unsigned int vD = 0; vD < LabelImageType::ImageDimension; vD++) {
                vLength += vOff[vD] * vOff[vD];
            }
            if(vLength > 0) {
                m_lengthWeightMask[vN] = 1.0f / sqrt(vLength);
            } else {
                m_lengthWeightMask[vN] = 0;
            }
        }
        
        m_NeedToUseBoundaryCondition = true;
    }
    
    template<typename TLabelImage, typename TEnergyDifference>
    void RCDiscreteContourLengthApproxEnergyEnergy<TLabelImage, TEnergyDifference>
    ::PrintSelf(std::ostream & os, Indent indent) const
    {
        Superclass::PrintSelf(os, indent);
    }
    
    template<typename TLabelImage, typename TEnergyDifference>
    inline 
    typename RCDiscreteContourLengthApproxEnergyEnergy<TLabelImage, TEnergyDifference>::
    InternalEnergyReturnType
    RCDiscreteContourLengthApproxEnergyEnergy<TLabelImage, TEnergyDifference>
    ::EvaluateAtIndexUsingLabel(IndexType const & aIndex, LabelPixelType aL) const {

        SizeType vSize;
        vSize.Fill(m_Radius);
        
        typedef ConstantBoundaryCondition<LabelImageType> ConstBoundaryCondType;
        ConstNeighborhoodIterator<LabelImageType, ConstBoundaryCondType> vImgIt(
                vSize, this->m_LabelImage, 
                this->m_LabelImage->GetBufferedRegion());
        vImgIt.SetNeedToUseBoundaryCondition(m_NeedToUseBoundaryCondition);
                
        vImgIt.SetLocation(aIndex);
        
        LabelPixelType vO; // abs pixel value around the center
        LabelPixelType vL = abs(aL); // abs pixel value in the center
        InternalEnergyReturnType vPenalty = 0;

        for (unsigned int vI = 0; vI < m_Offsets.size(); vI++) {

            vO = abs(vImgIt.GetPixel(m_Offsets[vI]));
            vPenalty += m_lengthWeightMask[vI] * H(abs(vL - vO));

        }
//        if(vPenalty == 2*vDim)
//            return itk::NumericTraits<float>::max();
        return vPenalty;
    }
        
    template<typename TLabelImage, typename TEnergyDifference>
    inline 
    typename RCDiscreteContourLengthApproxEnergyEnergy<TLabelImage, TEnergyDifference>::
    InternalEnergyReturnType
    RCDiscreteContourLengthApproxEnergyEnergy<TLabelImage, TEnergyDifference>
    ::EvaluateEnergyDifference(IndexType aInd, 
            LabelPixelType aLabelFrom, 
            LabelPixelType aLabelTo) {
       return this->m_Coefficient * (
               EvaluateAtIndexUsingLabel(aInd, aLabelTo) 
               - EvaluateAtIndexUsingLabel(aInd, aLabelFrom));
    }
    
    template<typename TLabelImage, typename TEnergyDifference>
    inline void 
    RCDiscreteContourLengthApproxEnergyEnergy<TLabelImage, TEnergyDifference>
    ::AddPoint(IndexType aInd, LabelAbsPixelType aLabel) {

    }
    
    template<typename TLabelImage, typename TEnergyDifference>
    inline void 
    RCDiscreteContourLengthApproxEnergyEnergy<TLabelImage, TEnergyDifference>
    ::RemovePoint(IndexType aInd, LabelAbsPixelType aLabel) {

    }
    
    template<typename TLabelImage, typename TEnergyDifference>
    void 
    RCDiscreteContourLengthApproxEnergyEnergy<TLabelImage, TEnergyDifference>
    ::KillRegion(LabelAbsPixelType aLabel) {

    }
        

} // end namespace itk
