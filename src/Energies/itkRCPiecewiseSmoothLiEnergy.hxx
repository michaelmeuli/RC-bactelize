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

#include "itkRCPiecewiseSmoothLiEnergy.h"

namespace itk
{
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    RCPiecewiseSmoothLiEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::RCPiecewiseSmoothLiEnergy() {
        m_Sigma = 2;
        m_GaussianImageSource = GaussianImageSourceType::New();
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void RCPiecewiseSmoothLiEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::PrintSelf(std::ostream & os, Indent indent) const
    {
        Superclass::PrintSelf(os, indent);
        os << indent << "Sigma: " << m_Sigma << std::endl;
    }
    
        
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void RCPiecewiseSmoothLiEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::PrepareEnergyCalculation() {
        
        typename GaussianImageSourceType::ArrayType vMean;
        typename GaussianImageSourceType::ArrayType vSigma;
        typename GaussianImageSourceType::SizeType vGaussImageSize;
        
        for(unsigned int vD = 0; vD < LabelImageType::ImageDimension;vD++) {
            float vScaledSigma = m_Sigma *
            this->m_DataImage->GetSpacing()[0] /
            this->m_DataImage->GetSpacing()[vD];
                        
            vMean[vD] = (vScaledSigma * 3);
            vSigma[vD] = (vScaledSigma);
            vGaussImageSize.Fill(vScaledSigma * 6 + 1);
        }
        
        m_GaussianImageSource->SetMean(vMean);
        m_GaussianImageSource->SetSigma(vSigma);
        m_GaussianImageSource->SetSize(vGaussImageSize);
        m_GaussianImageSource->SetNormalized(true);
        m_GaussianImageSource->Update();
        
        // normalization by ITK doesn't work??? center pixel has value 4.5 with
        // sigma = 3 ???
        ImageRegionIterator<DataImageType> vGaussIt(
        m_GaussianImageSource->GetOutput(), m_GaussianImageSource->GetOutput()->GetLargestPossibleRegion());
        DataPixelType vSum = 0;
        for(vGaussIt.GoToBegin(); !vGaussIt.IsAtEnd(); ++vGaussIt) {
            vSum += vGaussIt.Get();
        }
        for(vGaussIt.GoToBegin(); !vGaussIt.IsAtEnd(); ++vGaussIt) {
            vGaussIt.Set(vGaussIt.Get() / vSum);
        }
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    inline 
    typename RCPiecewiseSmoothLiEnergy<TLabelImage, TDataImage, TEnergyDifference>::
    ExternalEnergyReturnType
    RCPiecewiseSmoothLiEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::EvaluateEnergyDifference_(IndexType aInd, 
            LabelPixelType aLabelFrom, 
            LabelPixelType aLabelTo, 
            DataPixelType aImgValue) {
        
        typedef typename LabelImageType::RegionType RegionType;
        RegionType vRegion =
                m_GaussianImageSource->GetOutput()->GetLargestPossibleRegion();

        typedef typename LabelImageType::OffsetType OffsetType;
        OffsetType vOffset;
        for (unsigned int vD = 0; vD < LabelImageType::ImageDimension; vD++) {
            // flooring is on purpose!
            vOffset[vD] = (vRegion.GetSize())[vD] / 2;
        }


        vRegion.SetIndex(aInd - vOffset);
        // apply translation to the Gaussian filter
        m_GaussianImageSource->GetOutput()->SetBufferedRegion(vRegion);

        // After the cropping at the data-image boundaries, vRegion will be the
        // region to treat in the data-image space.
        vRegion.Crop(this->m_DataImage->GetBufferedRegion());
        vRegion.Crop(this->m_LabelImage->GetBufferedRegion()); // not necessary?

        ImageRegionConstIterator<LabelImageType> vLabelIt(this->m_LabelImage, vRegion);
        ImageRegionConstIterator<DataImageType> vDataIt(this->m_DataImage, vRegion);
        ImageRegionConstIteratorWithIndex<DataImageType> vGaussIt(
        m_GaussianImageSource->GetOutput(), vRegion);

        // we tentatively change the label value of the pixel of interest, such
        // that it will not be considered in the calculation.
//        LabelPixelType vStoredLabel = this->m_LabelImage->GetPixel(aInd);
//        this->m_LabelImage->SetPixel(aInd, aLabelFrom+aLabelTo);
        unsigned int vSubImgSize = 1;
        for(unsigned int vD = 0; vD < LabelImageType::ImageDimension; vD++){
            vSubImgSize *= vRegion.GetSize()[vD];
        }
        std::vector<LabelPixelType> vLabelSubImage;
        vLabelSubImage.resize(vSubImgSize);
        
        EnergyDifferenceType vWeightedSumFrom = 0;
        EnergyDifferenceType vWeightedSumTo = 0;
        EnergyDifferenceType vSumFrom = 0; 
        EnergyDifferenceType vSumTo = 0;
        EnergyDifferenceType vSquaredSumFrom = 0;
        EnergyDifferenceType vSquaredSumTo = 0;
        EnergyDifferenceType vNFrom = 0;
        EnergyDifferenceType vNTo = 0;
        EnergyDifferenceType vNormalizerFrom = 0;
        EnergyDifferenceType vNormalizerTo = 0;
        EnergyDifferenceType vGaussCenterWeight = 0;
        
        int vSubImgInd = -1;
        unsigned int vSubImgCenterIndex = 0;
        for (vLabelIt.GoToBegin(), vDataIt.GoToBegin(), vGaussIt.GoToBegin();
                !vLabelIt.IsAtEnd();
                ++vLabelIt, ++vDataIt, ++vGaussIt) {

            // copy the label image 
            vSubImgInd++;
            vLabelSubImage[vSubImgInd] = static_cast<unsigned int> (abs(vLabelIt.Get()));
            
            // figure out where the center pixel is in the (cropped) mask
            bool vCheckZero = true;
            for (unsigned int vD = 0; vD < LabelImageType::ImageDimension; vD++) {
                vCheckZero = vCheckZero && (aInd[vD] == vGaussIt.GetIndex()[vD]);
            }
            if(vCheckZero) {
                vGaussCenterWeight = vGaussIt.Get();
                vSubImgCenterIndex = vSubImgInd;
                continue;
            }
            
            // get weighted stats within mask; variance is not weighted
            if (static_cast<unsigned int> (abs(vLabelIt.Get())) == aLabelFrom) {
                vSumFrom += vDataIt.Get();
                vSquaredSumFrom += vDataIt.Get() * vDataIt.Get();
                double vV = vDataIt.Get() * vGaussIt.Get();
                vWeightedSumFrom += vV;
                vNormalizerFrom += vGaussIt.Get();
                vNFrom++;
            } else if (static_cast<unsigned int> (abs(vLabelIt.Get())) == aLabelTo) {
                vSumTo += vDataIt.Get();
                vSquaredSumTo += vDataIt.Get() * vDataIt.Get();
                double vV = vDataIt.Get() * vGaussIt.Get();
                vWeightedSumTo += vV;
                vNormalizerTo += vGaussIt.Get();
                vNTo++;
            }

        }

        EnergyDifferenceType vF1withX = (vWeightedSumFrom + aImgValue * vGaussCenterWeight) /
                (vNormalizerFrom + vGaussCenterWeight);
        EnergyDifferenceType vF1withoutX = (vNormalizerFrom == 0)?
            NumericTraits<EnergyDifferenceType>::max() : vWeightedSumFrom / vNormalizerFrom;
        
        EnergyDifferenceType vF2withX = (vWeightedSumTo + aImgValue * vGaussCenterWeight) /
                (vNormalizerTo + vGaussCenterWeight);
        EnergyDifferenceType vF2withoutX = (vNormalizerTo == 0)? 
            NumericTraits<EnergyDifferenceType>::max():vWeightedSumTo / vNormalizerTo;

        EnergyDifferenceType vE1_withX = 0;
        EnergyDifferenceType vE1_withoutX = 0;
        EnergyDifferenceType vE2_withX = 0;
        EnergyDifferenceType vE2_withoutX = 0;

        // Calculate the energy (convolution of the local CV model with a
        // Gaussian) considering the pixel of interest to be in 'from' label.
//        this->m_LabelImage->SetPixel(aInd, aLabelFrom);
        vLabelSubImage[vSubImgCenterIndex] = aLabelFrom;
        for (vSubImgInd = 0, vDataIt.GoToBegin(), vGaussIt.GoToBegin();
                !vDataIt.IsAtEnd();
                ++vSubImgInd, ++vDataIt, ++vGaussIt) {

            if (vLabelSubImage[vSubImgInd] == aLabelFrom) {
                vE1_withX += (vDataIt.Get() - vF1withX) *
                        (vDataIt.Get() - vF1withX) * vGaussIt.Get();
            } else if (vLabelSubImage[vSubImgInd] == aLabelTo) {
                vE2_withoutX += (vDataIt.Get() - vF2withoutX) *
                        (vDataIt.Get() - vF2withoutX) * vGaussIt.Get();
            }
        }
        // Calculate the energy considering the pixel of interest to have changed
        // to the new label 'to'.
//        this->m_LabelImage->SetPixel(aInd, aLabelTo);
        vLabelSubImage[vSubImgCenterIndex] = aLabelTo;
        for (vSubImgInd = 0, vDataIt.GoToBegin(), vGaussIt.GoToBegin();
                !vDataIt.IsAtEnd();
                ++vSubImgInd, ++vDataIt, ++vGaussIt) {

            if (vLabelSubImage[vSubImgInd] == aLabelFrom) {
                vE1_withoutX += (vDataIt.Get() - vF1withoutX) *
                        (vDataIt.Get() - vF1withoutX) * vGaussIt.Get();
            } else if (vLabelSubImage[vSubImgInd] == aLabelTo) {
                vE2_withX += (vDataIt.Get() - vF2withX) *
                        (vDataIt.Get() - vF2withX) * vGaussIt.Get();
            }
        }

        // set the value back to the original label; leave the labeling function
        // unchanged.
//        this->m_LabelImage->SetPixel(aInd, vStoredLabel);

        ExternalEnergyReturnType vRet;
        // return the energy gradient (the new state minus the current state)
        vRet.first = (vE1_withoutX + vE2_withX) - (vE1_withX + vE2_withoutX);
        vRet.first *= this->m_Coefficient;
        
        // local merging!
        if(vNFrom == 0 || vNTo == 0) {
            vRet.second = false;
        }else{
            vRet.second = this->CalculateKullbackLeiblerDistance(vSumFrom / vNFrom,
                    vSumTo / vNTo, 
                    this->CalculateVariance(vSquaredSumFrom, vSumFrom / vNFrom, vNFrom),
                    this->CalculateVariance(vSquaredSumTo / vNTo, vSumTo/vNTo, vNTo),
                    vNFrom, vNTo) 
                    < this->m_RegionMergingThreshold;
        }
        return vRet;
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    inline void 
    RCPiecewiseSmoothLiEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::AddPoint_(IndexType aInd, LabelAbsPixelType aLabel, DataPixelType aVal) {

    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    inline void 
    RCPiecewiseSmoothLiEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::RemovePoint_(IndexType aInd /* = 0 */, LabelAbsPixelType aLabel, DataPixelType aVal) {

    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void 
    RCPiecewiseSmoothLiEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::KillRegion_(LabelAbsPixelType aLabel) {

    }
        

} // end namespace itk