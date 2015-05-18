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

#include "itkRCConstantOutwardFlow.h"
#include "itkRCEnergyBaseClass.h"


namespace itk
{
    template<typename TLabelImage, typename TEnergyDifference>
    RCConstantOutwardFlowEnergy<TLabelImage, TEnergyDifference>
    ::RCConstantOutwardFlowEnergy() {
    }
    
    template<typename TLabelImage, typename TEnergyDifference>
    void RCConstantOutwardFlowEnergy<TLabelImage, TEnergyDifference>
    ::PrintSelf(std::ostream & os, Indent indent) const
    {
        Superclass::PrintSelf(os, indent);
    }
    
        
    template<typename TLabelImage, typename TEnergyDifference>
    inline 
    typename RCConstantOutwardFlowEnergy<TLabelImage, TEnergyDifference>::
    InternalEnergyReturnType
    RCConstantOutwardFlowEnergy<TLabelImage, TEnergyDifference>
    ::EvaluateEnergyDifference(IndexType aInd, 
            LabelPixelType aLabelFrom, 
            LabelPixelType aLabelTo) {
        
        InternalEnergyReturnType vEnergy = 0;
        if (aLabelFrom == 0) { // growing
            vEnergy -= this->m_Coefficient;
        } else if (aLabelTo == 0) { // shrinking
            vEnergy += this->m_Coefficient;
        }
       return (vEnergy);

    }
    
    template<typename TLabelImage, typename TEnergyDifference>
    inline void 
    RCConstantOutwardFlowEnergy<TLabelImage, TEnergyDifference>
    ::AddPoint(IndexType aInd, LabelAbsPixelType aLabel) {

    }
    
    template<typename TLabelImage, typename TEnergyDifference>
    inline void 
    RCConstantOutwardFlowEnergy<TLabelImage, TEnergyDifference>
    ::RemovePoint(IndexType aInd, LabelAbsPixelType aLabel) {

    }
    
    template<typename TLabelImage, typename TEnergyDifference>
    void 
    RCConstantOutwardFlowEnergy<TLabelImage, TEnergyDifference>
    ::KillRegion(LabelAbsPixelType aLabel) {

    }
        

} // end namespace itk