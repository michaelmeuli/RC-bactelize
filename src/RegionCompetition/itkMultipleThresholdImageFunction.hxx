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
#ifndef __itkMultipleThresholdImageFunction_hxx
#define __itkMultipleThresholdImageFunction_hxx

#include "itkMultipleThresholdImageFunction.h"

namespace itk {

    template <class TInputImage, class TCoordRep>
    MultipleThresholdImageFunction<TInputImage, TCoordRep>
    ::MultipleThresholdImageFunction() {
        m_NThresholds = 0;
    }

    /**
     * The values less than or equal to the value are inside
     */
    template <class TInputImage, class TCoordRep>
    void
    MultipleThresholdImageFunction<TInputImage, TCoordRep>
    ::AddThresholdBetween(PixelType lower, PixelType upper) {
        this->Modified();
        m_Thresholds.push_back(std::pair<PixelType, PixelType > (lower, upper));
        m_NThresholds += 1;
    }

    template <class TInputImage, class TCoordRep>
    void
    MultipleThresholdImageFunction<TInputImage, TCoordRep>
    ::PrintSelf(std::ostream& os, Indent indent) const {
        Superclass::PrintSelf(os, indent);
        os << indent << "Threshold: ";
        for (unsigned int vI = 0; vI < m_Thresholds.size(); vI++)
            os << vI + 1 << ": [" << m_Thresholds[vI].first << ", " << m_Thresholds[vI].second << "]";
        os << std::endl;
    }

} // end namespace itk

#endif
