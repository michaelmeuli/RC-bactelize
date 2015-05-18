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

#include "itkRCPiecewiseConstantSquareDistEnergy.h"

namespace itk
{
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    RCPiecewiseConstantSquareDistEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::RCPiecewiseConstantSquareDistEnergy() {

    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void RCPiecewiseConstantSquareDistEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::PrintSelf(std::ostream & os, Indent indent) const
    {
        Superclass::PrintSelf(os, indent);
        os << indent << "Sum statistic size: " << this->m_Sums.size() << std::endl;
        os << indent << "Count statistics size: " << this->m_Count.size() << std::endl;
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    inline 
    typename RCPiecewiseConstantSquareDistEnergy<TLabelImage, TDataImage, TEnergyDifference>::
    ExternalEnergyReturnType
    RCPiecewiseConstantSquareDistEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::EvaluateEnergyDifference_(IndexType aInd, 
            LabelPixelType aLabelFrom, 
            LabelPixelType aLabelTo, 
            DataPixelType aImgValue) {
        /*
         * Here we have the possibility to either put the current pixel
         * value to the BG, calculate the BG-mean and then calculate the
         * squared distance of the pixel to both means, BG and the mean
         * of the region (where the pixel currently still belongs to).
         *
         * The second option is to remove the pixel from the region and
         * calculate the new mean of this region. Then compare the squared
         * distance to both means. This option needs a region to be larger
         * than 1 pixel/voxel.
         * 
         */
        unsigned int vNTo = this->m_Count[aLabelTo];
        unsigned int vNFrom = this->m_Count[aLabelFrom];
        EnergyDifferenceType vNewToMean = (this->m_Sums[aLabelTo] + aImgValue) / (vNTo + 1);
        EnergyDifferenceType vMeansFrom = this->m_Sums[aLabelFrom] / vNFrom;
        EnergyDifferenceType vNewVarTo = this->CalculateVariance(
                this->m_Sums_2[aLabelTo] + aImgValue * aImgValue, vNewToMean, vNTo+1);
        EnergyDifferenceType vVarFrom = this->CalculateVariance(
                this->m_Sums_2[aLabelFrom], vMeansFrom, vNFrom);
        
        ExternalEnergyReturnType vRetType;
        vRetType.first = (aImgValue - vNewToMean) * (aImgValue - vNewToMean) -
                (aImgValue - vMeansFrom) * (aImgValue - vMeansFrom);
        vRetType.first *= Superclass::m_Coefficient;

        vRetType.second = this->CalculateKullbackLeiblerDistance(
                vNewToMean,vMeansFrom,vNewVarTo,vVarFrom,vNTo,vNFrom) < 
                Superclass::m_RegionMergingThreshold;
        
        return vRetType;
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    inline void 
    RCPiecewiseConstantSquareDistEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::AddPoint_(IndexType aInd, LabelAbsPixelType aLabel, DataPixelType aVal) {

    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    inline void 
    RCPiecewiseConstantSquareDistEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::RemovePoint_(IndexType aInd /* = 0 */, LabelAbsPixelType aLabel, DataPixelType aVal) {

    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void 
    RCPiecewiseConstantSquareDistEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::KillRegion_(LabelAbsPixelType aLabel) {

    }
        
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    typename RCPiecewiseConstantSquareDistEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::EnergyDifferenceType
    RCPiecewiseConstantSquareDistEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::CalculateTotalEnergy() {
        EnergyDifferenceType vEnergy = 0;
        
        typename DataImageType::Pointer vReconstructedImage;
        void * vimg = &vReconstructedImage;
        GenerateReconstructedImage(vimg);
        ImageRegionConstIteratorWithIndex<DataImageType> vDataIt(this->m_DataImage,
                this->m_DataImage->GetBufferedRegion());
        ImageRegionConstIterator<DataImageType> vReconstIt(vReconstructedImage,
                vReconstructedImage->GetBufferedRegion());
        for (vReconstIt.GoToBegin(), vDataIt.GoToBegin(); !vReconstIt.IsAtEnd();
                ++vReconstIt, ++vDataIt) {
            EnergyDifferenceType vDiff = static_cast<EnergyDifferenceType>(vReconstIt.Get()) - 
                    static_cast<EnergyDifferenceType>(vDataIt.Get());
            vEnergy += this->m_Coefficient * vDiff * vDiff;
        }
        return vEnergy;
    }
      
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void
    RCPiecewiseConstantSquareDistEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::GenerateReconstructedImage(void * aPointerToImg) {
        typename ReconstructedImageType::Pointer vReconstImage = ReconstructedImageType::New();

        vReconstImage->SetRequestedRegion(this->m_DataImage->GetRequestedRegion());
        vReconstImage->SetBufferedRegion(this->m_DataImage->GetBufferedRegion());
        vReconstImage->SetLargestPossibleRegion(this->m_DataImage->GetLargestPossibleRegion());
        vReconstImage->Allocate();
        vReconstImage->SetSpacing(this->m_DataImage->GetSpacing());

        /* reconstruct*/
        ImageRegionConstIteratorWithIndex<LabelImageType> vPixelIt(
                    this->m_LabelImage, this->m_LabelImage->GetBufferedRegion());
        
        for (vPixelIt.GoToBegin(); !vPixelIt.IsAtEnd(); ++vPixelIt) {
            LabelAbsPixelType vAbsLabel = abs(vPixelIt.Get());
            IndexType vCenterIndex = vPixelIt.GetIndex();
            vReconstImage->SetPixel(vCenterIndex, this->m_Sums[vAbsLabel] / this->m_Count[vAbsLabel]);
        }
        
        /* set the smart pointer */
        typename ReconstructedImageType::Pointer* vReturnPointer =
        (typename ReconstructedImageType::Pointer*) aPointerToImg;
        *vReturnPointer = vReconstImage;
        return;
    }
} // end namespace itk
