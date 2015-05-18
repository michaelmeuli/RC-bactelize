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
#include "itkRCParticleInteractionEnergyBaseClass.h"


namespace itk
{
    template <class TLabelImage, class TOptimizer, class TEnergyDifference>
    RCParticleInteractionEnergyBaseClass<TLabelImage, TOptimizer, TEnergyDifference>
    ::RCParticleInteractionEnergyBaseClass() {//: RCEnergyBaseClass() {
        m_NeighborsOffsets_BG_Connectivity =
                BackgroundConnectivityType::GetInstance().GetNeighborsITKOffsets();
        m_NeighborhoodSize_BG_Connectivity =
                BackgroundConnectivityType::GetInstance().GetNumberOfNeighbors();
        m_ParticleInteractionRadius = 1;
    }
    
    template <class TLabelImage, class TOptimizer, class TEnergyDifference>
    void RCParticleInteractionEnergyBaseClass<TLabelImage, TOptimizer, TEnergyDifference>
    ::PrintSelf(std::ostream & os, Indent indent) const
    {
        Superclass::PrintSelf(os, indent);
        os << indent << "ParticleInteractionRadius: " << m_ParticleInteractionRadius << std::endl;
    }
    
    template <class TLabelImage, class TOptimizer, class TEnergyDifference>
    typename RCParticleInteractionEnergyBaseClass<TLabelImage, TOptimizer, TEnergyDifference>
    ::EnergyDifferenceType 
    RCParticleInteractionEnergyBaseClass<TLabelImage, TOptimizer, TEnergyDifference>
    ::ApproximateSurfaceNormal(IndexType aIndex, VectorType& aNormal) {

        for (unsigned int vD = 0; vD < LabelImageType::ImageDimension; vD++) {
            aNormal[vD] = 0;
        }

        LabelAbsPixelType vLabelOfCenter = abs(this->m_LabelImage->GetPixel(aIndex));
        for (unsigned int vN = 0; vN < m_NeighborhoodSize_BG_Connectivity; vN++) {
            OffsetType vOff = m_NeighborsOffsets_BG_Connectivity[vN];
            if(abs(this->m_LabelImage->GetPixel(aIndex + vOff)) != vLabelOfCenter) {
                for (unsigned int vD = 0; vD < LabelImageType::ImageDimension; vD++) {
                    aNormal[vD] += vOff[vD];
                }
            }
        }
        
        /// Normalize the vector
        float vNorm = 0;
        for (unsigned int vD = 0; vD < LabelImageType::ImageDimension; vD++) {
            vNorm += aNormal[vD] * aNormal[vD];
        }
        vNorm = sqrt(vNorm);
        
        return vNorm;
    }
    

        
//    template <class TLabelImage, class TOptimizer, class TEnergyDifference>
//    void RCParticleInteractionEnergyBaseClass<TLabelImage, TOptimizer, TEnergyDifference>
//    ::EvaluateEnergyDifferences(const std::vector<IndexType>& aInds, 
//            const std::vector<LabelPixelType>& aLabelBefore, 
//            const std::vector<LabelPixelType>& aLabelAfter,
//            std::vector<ParticleInteractionEnergyReturnType>& aReturnVec) {
//        
//        aReturnVec.clear();
//        aReturnVec.resize(aInds.size());
//        
//        for(unsigned int vI = 0; vI < aReturnVec.size(); vI++) {
//            aReturnVec[vI] =
//                    EvaluateEnergyDifference(aInds[vI], 
//                    aLabelBefore[vI], 
//                    aLabelAfter[vI]);
//        }
//    }
    
} // end namespace itk
