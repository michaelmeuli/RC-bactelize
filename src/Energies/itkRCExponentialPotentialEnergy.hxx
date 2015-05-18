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

#include "itkRCExponentialPotentialEnergy.h"


namespace itk
{
    template<typename TLabelImage, typename TOptimizer, typename TEnergyDifference>
    RCExponentialPotentialEnergy<TLabelImage, TOptimizer, TEnergyDifference>
    ::RCExponentialPotentialEnergy() {

    }
    
    template<typename TLabelImage, typename TOptimizer, typename TEnergyDifference>
    void RCExponentialPotentialEnergy<TLabelImage, TOptimizer, TEnergyDifference>
    ::PrintSelf(std::ostream & os, Indent indent) const
    {
        Superclass::PrintSelf(os, indent);
    }
    
    
        
    template<typename TLabelImage, typename TOptimizer, typename TEnergyDifference>
    inline 
    typename RCExponentialPotentialEnergy<TLabelImage, TOptimizer, TEnergyDifference>::
    ParticleInteractionEnergyReturnType
    RCExponentialPotentialEnergy<TLabelImage, TOptimizer, TEnergyDifference>
    ::EvaluateEnergyDifference(IndexType aInd, 
            InteractionParticlesIteratorType aBegin, InteractionParticlesIteratorType aEnd) {

        ParticleInteractionEnergyReturnType vGradientFlow = 0;
        
        VectorType vNormal;
        float vNorm = this->ApproximateSurfaceNormal(aInd, vNormal);
        LabelAbsPixelType vLabelOfCenter = static_cast<LabelAbsPixelType>(abs(this->m_LabelImage->GetPixel(aInd)));
                
        if (vNorm > 0.001f) {
            for (; aBegin != aEnd; ++aBegin) {

                if (aBegin->second.m_label != vLabelOfCenter) {
                    continue;
                }

                // euclediean distance between p and q
                float vSqDist = 0;
                for (unsigned int vD = 0; vD < LabelImageType::ImageDimension; vD++) {
                    vSqDist += (aBegin->first[vD] - aInd[vD]) * (aBegin->first[vD] - aInd[vD]);
                }

                /**
                 * Exponential potential
                 * In this case, the interaction radius is the standard dev. of
                 * the corresponding Gaussian.
                 */
                // MATLAB: -2./(sigma^3*2*pi).*(p(i,:)-p(j,:)).*exp(-(d).^2/(sigma.^2));
                // MATLAB: dvE_dx = dE_dx*v(i,:)'
                // we leave away the 2/sqrt(2* pi) term since it is constant
                float vPISigma2 = this->m_ParticleInteractionRadius * this->m_ParticleInteractionRadius;
                float vExpTerm = -1.0f /
                        (vPISigma2) *
                        exp(-vSqDist / vPISigma2);

                for (unsigned int vD = 0; vD < LabelImageType::ImageDimension; vD++) {
                    vGradientFlow += vExpTerm * (aInd[vD] - aBegin->first[vD]) * vNormal[vD];
                }
            }
        }
        return this->m_Coefficient * (vGradientFlow);
    }
    
    template <class TLabelImage, class TOptimizer, class TEnergyDifference>
    typename RCExponentialPotentialEnergy<TLabelImage, TOptimizer, TEnergyDifference>
    ::SizeType 
    RCExponentialPotentialEnergy<TLabelImage, TOptimizer, TEnergyDifference>
    ::GetMinimalIteratorSize(){
        SizeType vSize;
        vSize.Fill(3 * this->m_ParticleInteractionRadius);
        return vSize;
    }
    
    template<typename TLabelImage, typename TOptimizer, typename TEnergyDifference>
    inline void 
    RCExponentialPotentialEnergy<TLabelImage, TOptimizer, TEnergyDifference>
    ::AddPoint(IndexType aInd, LabelAbsPixelType aLabel) {

    }
    
    template<typename TLabelImage, typename TOptimizer, typename TEnergyDifference>
    inline void 
    RCExponentialPotentialEnergy<TLabelImage, TOptimizer, TEnergyDifference>
    ::RemovePoint(IndexType aInd, LabelAbsPixelType aLabel) {

    }
    
    template<typename TLabelImage, typename TOptimizer, typename TEnergyDifference>
    void 
    RCExponentialPotentialEnergy<TLabelImage, TOptimizer, TEnergyDifference>
    ::KillRegion(LabelAbsPixelType aLabel) {

    }
        

} // end namespace itk