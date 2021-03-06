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

#include "itkRCCenterOfMassDirectedFlow.h"


namespace itk
{
    template<typename TLabelImage, typename TOptimizer, typename TEnergyDifference>
    RCCenterOfMassDirectedFlow<TLabelImage, TOptimizer, TEnergyDifference>
    ::RCCenterOfMassDirectedFlow() {

    }
    
    template<typename TLabelImage, typename TOptimizer, typename TEnergyDifference>
    void RCCenterOfMassDirectedFlow<TLabelImage, TOptimizer, TEnergyDifference>
    ::PrintSelf(std::ostream & os, Indent indent) const
    {
        Superclass::PrintSelf(os, indent);
    }
    
    
        
    template<typename TLabelImage, typename TOptimizer, typename TEnergyDifference>
    inline 
    typename RCCenterOfMassDirectedFlow<TLabelImage, TOptimizer, TEnergyDifference>::
    ParticleInteractionEnergyReturnType
    RCCenterOfMassDirectedFlow<TLabelImage, TOptimizer, TEnergyDifference>
    ::EvaluateEnergyDifference(IndexType aInd, 
            InteractionParticlesIteratorType aBegin, InteractionParticlesIteratorType aEnd) {

        ParticleInteractionEnergyReturnType vFlow = 0;
        float vSpringConstant_2 = this->m_ParticleInteractionRadius * this->m_ParticleInteractionRadius;
        
        VectorType vNormal;
        float vNorm = this->ApproximateSurfaceNormal(aInd, vNormal);
        LabelAbsPixelType vLabelOfCenter = static_cast<LabelAbsPixelType>(abs(this->m_LabelImage->GetPixel(aInd)));
                
        VectorType vCenterOfMass;
        vCenterOfMass.Fill(0);
        
        
        if (vNorm > 0.001f) {
            unsigned int vNeigh = 0;
            for (; aBegin != aEnd; ++aBegin) {

                if (aBegin->second.m_label != vLabelOfCenter) {
                    continue;
                }

                // euclediean distance between p and q
                float vSqDist = 0;
                for (unsigned int vD = 0; vD < LabelImageType::ImageDimension; vD++) {
                    vSqDist += (aBegin->first[vD] - aInd[vD]) * (aBegin->first[vD] - aInd[vD]);
                }
                
                /*
                 * Calculate the center of mass and project to normal
                 */
                if (vSqDist < this->m_ParticleInteractionRadius * this->m_ParticleInteractionRadius) {
                    vNeigh++;
                    for (unsigned int vD = 0; vD < LabelImageType::ImageDimension; vD++) {
                        vCenterOfMass[vD] = vCenterOfMass[vD] + aBegin->first[vD];
                    }
                }
                if (vNeigh > 0) {
                    for (unsigned int vD = 0; vD < LabelImageType::ImageDimension; vD++) {
                        vFlow += -(aInd[vD] - vCenterOfMass[vD] / vNeigh) * vNormal[vD];
                    }
                }
            }
        }
        return this->m_Coefficient * (vFlow);
    }
    
    template <class TLabelImage, class TOptimizer, class TEnergyDifference>
    typename RCCenterOfMassDirectedFlow<TLabelImage, TOptimizer, TEnergyDifference>
    ::SizeType 
    RCCenterOfMassDirectedFlow<TLabelImage, TOptimizer, TEnergyDifference>
    ::GetMinimalIteratorSize(){
        SizeType vSize;
        vSize.Fill(this->m_ParticleInteractionRadius);
        return vSize;
    }
        
    template<typename TLabelImage, typename TOptimizer, typename TEnergyDifference>
    inline void 
    RCCenterOfMassDirectedFlow<TLabelImage, TOptimizer, TEnergyDifference>
    ::AddPoint(IndexType aInd, LabelAbsPixelType aLabel) {

    }
    
    template<typename TLabelImage, typename TOptimizer, typename TEnergyDifference>
    inline void 
    RCCenterOfMassDirectedFlow<TLabelImage, TOptimizer, TEnergyDifference>
    ::RemovePoint(IndexType aInd, LabelAbsPixelType aLabel) {

    }
    
    template<typename TLabelImage, typename TOptimizer, typename TEnergyDifference>
    void 
    RCCenterOfMassDirectedFlow<TLabelImage, TOptimizer, TEnergyDifference>
    ::KillRegion(LabelAbsPixelType aLabel) {

    }
        

} // end namespace itk