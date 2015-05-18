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

#include "itkRCKybicKratkyContourLengthApproxEnergy.h"
#include "itkSphereBitmapImageSource.h"
#include "itkRCEnergyBaseClass.h"


namespace itk
{
    template<typename TLabelImage, typename TEnergyDifference>
    RCKybicKratkyContourLengthApproxEnergyEnergy<TLabelImage, TEnergyDifference>
    ::RCKybicKratkyContourLengthApproxEnergyEnergy() {
        // the following 2 members are correctly initialized in 
        // PrepareEnergyCalculation
        m_Prefactor = 0;
        m_Volume = 0;
        m_Radius = 4;
        m_SphereMaskForCurvature = SphereMaskImageSourceType::New();
    }
    
    template<typename TLabelImage, typename TEnergyDifference>
    void RCKybicKratkyContourLengthApproxEnergyEnergy<TLabelImage, TEnergyDifference>
    ::PrintSelf(std::ostream & os, Indent indent) const
    {
        Superclass::PrintSelf(os, indent);
        os << indent << "Radius: " << m_Radius << std::endl;
        //os << indent << "Prefactor: " << m_Prefactor << std::endl;
    }
    
    template<typename TLabelImage, typename TEnergyDifference>
    void RCKybicKratkyContourLengthApproxEnergyEnergy<TLabelImage, TEnergyDifference>
    ::PrepareEnergyCalculation() {
        
        m_Volume = this->CalculateScaledSphereVolume(m_Radius);
        if(LabelImageType::ImageDimension == 2) {
            m_Prefactor = 3.0 * 3.141592 / m_Radius;
        } else if(LabelImageType::ImageDimension == 3) {
            m_Prefactor = 16.0f / (3.0f * m_Radius);
        } else {
            itkExceptionMacro("This energy is only defined for 2 and 3 dimensional"
            "images.");
        }
            
        typename SphereMaskImageSourceType::SizeType vMaskSize;
        typename SphereMaskImageSourceType::ArrayType vRadius;
        for (unsigned int vD = 0; vD < LabelImageType::ImageDimension; vD++) {

            // the radius is expected to be given in px size of the first
            // axis. We scale for the all the following dimensions according
            // to the image spacing.
            float vScaledRadius = m_Radius *
            this->m_LabelImage->GetSpacing()[0] /
            this->m_LabelImage->GetSpacing()[vD];
            
            vRadius[vD] = vScaledRadius;
            vMaskSize[vD] = ceil(vScaledRadius) * 2 + 1;
        }
        m_SphereMaskForCurvature->SetRadius(vRadius);
        m_SphereMaskForCurvature->SetSize(vMaskSize);
        m_SphereMaskForCurvature->Update();
    }
        
    template<typename TLabelImage, typename TEnergyDifference>
    inline 
    typename RCKybicKratkyContourLengthApproxEnergyEnergy<TLabelImage, TEnergyDifference>::
    InternalEnergyReturnType
    RCKybicKratkyContourLengthApproxEnergyEnergy<TLabelImage, TEnergyDifference>
    ::EvaluateEnergyDifference(IndexType aInd, 
            LabelPixelType aLabelFrom, 
            LabelPixelType aLabelTo) {
        
        InternalEnergyReturnType vCurvatureFlow = 0;

        // vRegion is the size of our temporary window
        RegionType vRegion;
        // vOffset is basically the difference of the center and the start of the window
        OffsetType vOffset;

        //read out the size of the mask
        vRegion = m_SphereMaskForCurvature->GetOutput()->GetBufferedRegion();
        for (unsigned int vD = 0; vD < LabelImageType::ImageDimension; vD++) {
            // flooring is on purpose!
            vOffset[vD] = (vRegion.GetSize())[vD] / 2;
        }
        // translate the region
        vRegion.SetIndex(aInd - vOffset);
        // apply translation to the mask
        m_SphereMaskForCurvature->GetOutput()->SetBufferedRegion(vRegion);

        // we assume the label image to be in memory and to be as large as the
        // data image. 
        vRegion.Crop(this->m_LabelImage->GetBufferedRegion());

        ImageRegionConstIterator<LabelImageType> vLabelIt(this->m_LabelImage, vRegion);
        ImageRegionConstIterator<MaskImageType> vMaskIt(
                m_SphereMaskForCurvature->GetOutput(),
                vRegion);

        if (aLabelFrom == 0) { // && aTo != 0) { //growing
            unsigned int vN = 0;
            // This is a point on the contour (innerlist) OR
            // touching the contour (Outer list)
            for (vLabelIt.GoToBegin(), vMaskIt.GoToBegin();
                    !vLabelIt.IsAtEnd();
                    ++vLabelIt, ++vMaskIt) {
                if (vMaskIt.Get()) {
                    if (static_cast<unsigned int> (abs(vLabelIt.Get())) == aLabelTo) {
                        vN++;
                    }
                }
            }
            vCurvatureFlow -= m_Prefactor * (static_cast<float> (vN) / m_Volume - 0.5f);
        } else {
            if (aLabelTo == 0) {//proper shrinking
                unsigned int vN = 0;
                // This is a point on the contour (innerlist) OR
                // touching the contour (Outer list)
                for (vLabelIt.GoToBegin(), vMaskIt.GoToBegin();
                        !vLabelIt.IsAtEnd();
                        ++vLabelIt, ++vMaskIt) {
                    if (vMaskIt.Get()) {
                        if (static_cast<unsigned int> (abs(vLabelIt.Get())) == aLabelFrom) {
                            vN++;
                        }
                    }
                }
                vCurvatureFlow += m_Prefactor * (static_cast<float> (vN) / m_Volume - 0.5f);
                
            } else { // fighting fronts
                unsigned int vNFrom = 0;
                unsigned int vNTo = 0;
                for (vLabelIt.GoToBegin(), vMaskIt.GoToBegin();
                        !vLabelIt.IsAtEnd();
                        ++vLabelIt, ++vMaskIt) {
                    if (vMaskIt.Get()) {
                        if (static_cast<unsigned int> (abs(vLabelIt.Get())) == aLabelTo) {
                            vNTo++;
                        }
                        if (static_cast<unsigned int> (abs(vLabelIt.Get())) == aLabelFrom) {
                            vNFrom++;
                        }
                    }
                }
                vCurvatureFlow -= m_Prefactor * (static_cast<float> (vNTo) / m_Volume - 0.5f);
                vCurvatureFlow += m_Prefactor * (static_cast<float> (vNFrom) / m_Volume - 0.5f);
            }
        }

       return this->m_Coefficient * (vCurvatureFlow);
    }
    
    template<typename TLabelImage, typename TEnergyDifference>
    inline void 
    RCKybicKratkyContourLengthApproxEnergyEnergy<TLabelImage, TEnergyDifference>
    ::AddPoint(IndexType aInd, LabelAbsPixelType aLabel) {
    }
    
    template<typename TLabelImage, typename TEnergyDifference>
    inline void 
    RCKybicKratkyContourLengthApproxEnergyEnergy<TLabelImage, TEnergyDifference>
    ::RemovePoint(IndexType aInd, LabelAbsPixelType aLabel) {
    }
    
    template<typename TLabelImage, typename TEnergyDifference>
    void 
    RCKybicKratkyContourLengthApproxEnergyEnergy<TLabelImage, TEnergyDifference>
    ::KillRegion(LabelAbsPixelType aLabel) {

    }
    
//    template<typename TLabelImage, typename TEnergyDifference>
//    typename RCKybicKratkyContourLengthApproxEnergyEnergy<TLabelImage, TEnergyDifference>::
//    InternalEnergyReturnType 
//    RCKybicKratkyContourLengthApproxEnergyEnergy<TLabelImage, TEnergyDifference>
//    ::CalculateScaledSphereVolume(float aRadiusX) {
//
//        InternalEnergyReturnType vVolume;
//
//        float aRadiusY = aRadiusX *
//                this->m_LabelImage->GetSpacing()[0]
//                / this->m_LabelImage->GetSpacing()[1];
//
//        if (LabelImageType::ImageDimension == 2) {
//            vVolume = 3.141592f * aRadiusX * aRadiusY;
//        } else if (LabelImageType::ImageDimension == 3) {
//            float aRadiusZ = aRadiusX *
//                    this->m_LabelImage->GetSpacing()[0]
//                    / this->m_LabelImage->GetSpacing()[2];
//
//            vVolume = 4.0f / 3.0f * 3.141592f *
//                    aRadiusX * aRadiusY * aRadiusZ;
//        }
//
//        return vVolume;
//    }
        
} // end namespace itk
