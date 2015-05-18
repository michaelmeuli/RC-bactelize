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

#ifndef ITKRCUTILS_HXX
#define	ITKRCUTILS_HXX

#include "itkRCUtils.h"

namespace itk {
    
    template <class TLabelImage, unsigned int VCellDim>
    unsigned int RCUtils<TLabelImage, VCellDim>::m_NeighborhoodSize_FG_Connectivity =
    RCUtils<TLabelImage, VCellDim>::ForegroundConnectivityType::GetInstance().GetNumberOfNeighbors();
    
    template <class TLabelImage, unsigned int VCellDim>
    const typename RCUtils<TLabelImage, VCellDim>::OffsetType * RCUtils<TLabelImage, VCellDim>
    ::m_NeighborsOffsets_FG_Connectivity = 
    RCUtils<TLabelImage, VCellDim>::ForegroundConnectivityType::GetInstance().GetNeighborsITKOffsets();
    
    template <class TLabelImage, unsigned int VCellDim>
    unsigned int RCUtils<TLabelImage, VCellDim>::m_NeighborhoodSize_BG_Connectivity =
    RCUtils<TLabelImage, VCellDim>::BackgroundConnectivityType::GetInstance().GetNumberOfNeighbors();
    
    template <class TLabelImage, unsigned int VCellDim>
    const typename RCUtils<TLabelImage, VCellDim>::OffsetType * RCUtils<TLabelImage, VCellDim>
    ::m_NeighborsOffsets_BG_Connectivity =
    RCUtils<TLabelImage, VCellDim>::BackgroundConnectivityType::GetInstance().GetNeighborsITKOffsets();
    
    template <class TLabelImage, unsigned int VCellDim>
    inline bool RCUtils<TLabelImage, VCellDim>::
    IsEnclosedByLabel_BGConnectivity(LabelImagePointerType aLabelImage,
            const LabelImageIndexType& aIndex, 
            LabelAbsPixelType aAbsLabel) {
        
        for (unsigned int vI = 0; vI < m_NeighborhoodSize_BG_Connectivity; vI++) {
            if (static_cast<unsigned int> (abs(aLabelImage->GetPixel(
                    aIndex + m_NeighborsOffsets_BG_Connectivity[vI]))) != aAbsLabel)
                return false;
        }
        return true;
    }
    
    /**
     * The method checks if any of the FG-connected neighbors has the label aAbsLabel
     * @return true if no neighbor has the same label
     */
    template <class TLabelImage, unsigned int VCellDim>
    inline bool RCUtils<TLabelImage, VCellDim>::
    IsSingleFGPoint(LabelImagePointerType aLabelImage, 
            const LabelImageIndexType& aIndex, LabelAbsPixelType aAbsLabel) {

        for (unsigned int vI = 0; vI < m_NeighborhoodSize_FG_Connectivity; vI++) {
            if (static_cast<unsigned int> (abs(aLabelImage->GetPixel(
                    aIndex + m_NeighborsOffsets_FG_Connectivity[vI]))) == aAbsLabel) {
                return false;
            }
        }
        return true;
    }
    
    template <class TLabelImage, unsigned int VCellDim>
    inline bool RCUtils<TLabelImage, VCellDim>::
    IsBoundaryPoint(LabelImagePointerType aLabelImage, const LabelImageIndexType& aIndex)  {
        LabelAbsPixelType vLabelAbs = static_cast<unsigned int> (abs(aLabelImage->GetPixel(aIndex)));
        for (unsigned int vI = 0; vI < m_NeighborhoodSize_FG_Connectivity; vI++) {
            if (static_cast<unsigned int> (
                    abs(aLabelImage->GetPixel(aIndex + m_NeighborsOffsets_FG_Connectivity[vI])))
                    != vLabelAbs) {
                return true;
            }
        }
        return false;
    }
    
}

#endif