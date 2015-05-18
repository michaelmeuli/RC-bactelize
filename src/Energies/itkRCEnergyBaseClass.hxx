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
#include "itkRCEnergyBaseClass.h"

namespace itk
{
    template <class TLabelImage, class TEnergyDifference>
    RCEnergyBaseClass<TLabelImage, TEnergyDifference>
    ::RCEnergyBaseClass() {
        m_LabelImage = NULL;
    }
    
    template <class TLabelImage, class TEnergyDifference>
    void
    RCEnergyBaseClass<TLabelImage, TEnergyDifference>
    ::PrintSelf(std::ostream & os, Indent indent) const
    {
        Superclass::PrintSelf(os, indent);
        os << indent << "Energy term coefficient: " << m_Coefficient << std::endl;
//        os << indent << "Label image ptr: " << m_LabelImage << std::endl;
    }
    
    template <class TLabelImage, class TEnergyDifference>
    typename RCEnergyBaseClass<TLabelImage, TEnergyDifference>::EnergyDifferenceType 
    RCEnergyBaseClass<TLabelImage, TEnergyDifference>
    ::CalculateVariance(EnergyDifferenceType aSumSq, EnergyDifferenceType aMean, EnergyDifferenceType aN) {
        if(aN < 2) return 1; //TODO: what would be appropriate?
        return (aSumSq - aN * aMean * aMean)/(aN - 1);
    } 
    
    template <class TLabelImage, class TEnergyDifference>
    typename RCEnergyBaseClass<TLabelImage, TEnergyDifference>::EnergyDifferenceType 
    RCEnergyBaseClass<TLabelImage, TEnergyDifference>
    ::CalculateScaledSphereVolume(float aRadiusX) {

        EnergyDifferenceType vVolume;

        float aRadiusY = aRadiusX *
                m_LabelImage->GetSpacing()[0]
                / m_LabelImage->GetSpacing()[1];

        if (LabelImageType::ImageDimension == 2) {
            vVolume = 3.141592f * aRadiusX * aRadiusY;
        } else if (LabelImageType::ImageDimension == 3) {
            float aRadiusZ = aRadiusX *
                    m_LabelImage->GetSpacing()[0]
                    / m_LabelImage->GetSpacing()[2];

            vVolume = 4.0f / 3.0f * 3.141592f *
                    aRadiusX * aRadiusY * aRadiusZ;
        }

        return vVolume;
    }
} // end namespace itk
