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

#include "itkRCPiecewiseConstantDeconvolutionEnergyBaseClass.h"

namespace itk
{
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    RCPiecewiseConstantDeconvolutionEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
    ::RCPiecewiseConstantDeconvolutionEnergyBaseClass() {
        m_DeconvolutionModelImage = InternalImageType::New();
        m_PSF = InternalImageType::New();
        m_DoNormalizePSF = false;
        m_PSFInput = NULL;
        m_SwitchPointMode = false;
        m_OptimizationMode = false;
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void RCPiecewiseConstantDeconvolutionEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
    ::PrintSelf(std::ostream & os, Indent indent) const
    {
        Superclass::PrintSelf(os, indent);
        os << indent << "Do normalize PSF: " << m_DoNormalizePSF << std::endl;
        os << indent << "Optimization mode: " << m_OptimizationMode << std::endl;
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    inline void 
    RCPiecewiseConstantDeconvolutionEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
    ::AddPoint_(IndexType aInd, LabelAbsPixelType aLabel, DataPixelType aVal) {
        if(m_Intensities.find(aLabel) == m_Intensities.end()) {
            m_NewRegions.insert(aLabel);
        }
        if(!m_SwitchPointMode) { // in that case the label image will be updated directly
            // extend the standard statistics with the model image J:
            UpdateModelImage(aInd, 0, aLabel);
        }
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    inline void 
    RCPiecewiseConstantDeconvolutionEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
    ::RemovePoint_(IndexType aInd /* = 0 */, LabelAbsPixelType aLabel, DataPixelType aVal) {
        if(!m_SwitchPointMode) { // in that case the label image will be updated directly
            // extend the standard statistics with the model image J:
            UpdateModelImage(aInd, aLabel, 0);
        }
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    inline void 
    RCPiecewiseConstantDeconvolutionEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
    ::SwitchPoint(IndexType aInd, LabelAbsPixelType aLabelFrom, 
            LabelAbsPixelType aLabelTo, DataPixelType aVal){
        
        m_SwitchPointMode = true;
        Superclass::RemovePoint(aInd, aLabelFrom, aVal); // update standard stats
        Superclass::AddPoint(aInd, aLabelTo, aVal); // update standard stats
        UpdateModelImage(aInd, aLabelFrom, aLabelTo); 
        m_SwitchPointMode = false;
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void 
    RCPiecewiseConstantDeconvolutionEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
    ::KillRegion_(LabelAbsPixelType aLabel) {
        this->m_Intensities.erase(aLabel);
    }
    
        
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void 
    RCPiecewiseConstantDeconvolutionEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
    ::PrepareEnergyCalculation() {
        
        /* 
         * Prepare the PSF (check inputs and normalize)
         */
        if(!m_PSFInput) itkExceptionMacro("No PSF image provided.");
//        m_PSFInput->Update();
        
        EnergyDifferenceType vScale = 1;
        if(m_DoNormalizePSF) {
            typedef itk::StatisticsImageFilter<InternalImageType> StatisticsFilterType;
            typename StatisticsFilterType::Pointer vStatisticsFilter = StatisticsFilterType::New();
            vStatisticsFilter->SetInput(m_PSFInput);
            vStatisticsFilter->Update();
            vScale = vStatisticsFilter->GetSum();
        }
        
        // Use a shift-scale filter and use its output as our psf. This is safe
        // as the smart pointer will increase the reference count of the 
        // output image.
        typedef ShiftScaleImageFilter<DataImageType, InternalImageType> 
        ShiftScaleFilterType;
        typename ShiftScaleFilterType::Pointer vSSFilter = ShiftScaleFilterType::New();
        vSSFilter->SetInput(m_PSFInput);
        vSSFilter->SetShift(0);
        if(vScale == 0) itkExceptionMacro("PSF is empty.");
        vSSFilter->SetScale(static_cast<EnergyDifferenceType>(1.0/vScale));
        vSSFilter->Update();
        m_PSF = vSSFilter->GetOutput();
        
        // Initialize the intensities to the mean (as a first guess)
        SumStatisticsIteratorType vStatIt = this->m_Sums.begin();
        for(; vStatIt != this->m_Sums.end(); ++vStatIt) {
            LabelAbsPixelType vL = vStatIt->first;
            m_Intensities[vL] = vStatIt->second / this->m_Count[vL];
        }
        
        /// First, generate a rough estimate of the model image using the means
        /// as an intensity estimate. In a second step refine the estimates
        /// in RenewDeconvolutionStatistics().
        void * vimg = &m_DeconvolutionModelImage;
        GenerateModelImage(vimg, this->m_LabelImage, this->m_DataImage, m_PSF, &m_Intensities);
        RenewDeconvolutionStatistics(this->m_LabelImage, this->m_DataImage);
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void 
    RCPiecewiseConstantDeconvolutionEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
    ::PrepareEnergyCalculationForIteration() {
        RenewDeconvolutionStatistics(this->m_LabelImage, this->m_DataImage);
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void 
    RCPiecewiseConstantDeconvolutionEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
    ::RenewDeconvolutionStatistics(
    LabelImageConstPointerType aInitImage,
    DataImageConstPointerType aDataImage) {

        // In case of a new region (due to topological changes), we have to re-
        // estimate the model image before we can fit to the data.
        bool vAnyNewRegion = false;
        for(LabelAbsSetIteratorType vNewLabelsIt = m_NewRegions.begin(); 
                vNewLabelsIt != m_NewRegions.end(); ++vNewLabelsIt){
            
            this->m_Intensities[*vNewLabelsIt] = 
                    this->m_Sums[*vNewLabelsIt] / this->m_Count[*vNewLabelsIt];
            
            vAnyNewRegion = true;
        }
        if(vAnyNewRegion) {
            void * vimg = &m_DeconvolutionModelImage;
            GenerateModelImage(vimg, this->m_LabelImage, this->m_DataImage, m_PSF, &m_Intensities);
            m_NewRegions.clear();
        }
        /*
         * Generate a model image using rough estimates of the intensities. Here,
         * we use the old intensity values.
         * For all FG regions, find the median of the scaling factor.
         */

        // The BG region is not fitted above(since it may be very large and thus
        // the mean is a good approx), set it to the mean value:
        InternalPixelType vOldBG = m_Intensities[0];
        //        m_Intensities[0] = m_Means[0];

        typedef itk::SubtractImageFilter<InternalImageType, InternalImageType, InternalImageType>
                SubtractModelImageFilterType;
        typename SubtractModelImageFilterType::Pointer vModelImageSubtractionFilter =
                SubtractModelImageFilterType::New();
        vModelImageSubtractionFilter->SetInput1(m_DeconvolutionModelImage);
        vModelImageSubtractionFilter->SetConstant2(vOldBG); // here it has to be the old BG value!

        typedef itk::SubtractImageFilter<InternalImageType, InternalImageType, InternalImageType>
                SubtractInputImageFilterType;
        typename SubtractInputImageFilterType::Pointer vInputImageSubtractionFilter =
                SubtractInputImageFilterType::New();
        vInputImageSubtractionFilter->SetInput1(this->m_DataImage);
        vInputImageSubtractionFilter->SetConstant2(vOldBG);

        typedef itk::DivideImageFilter<InternalImageType, InternalImageType, InternalImageType>
                DivisionFilterType;
        typename DivisionFilterType::Pointer vDivisionFilter = DivisionFilterType::New();
        vDivisionFilter->SetInput1(vInputImageSubtractionFilter->GetOutput());
        vDivisionFilter->SetInput2(vModelImageSubtractionFilter->GetOutput());
        vDivisionFilter->Update();

        typedef itk::SubtractImageFilter<InternalImageType, InternalImageType, InternalImageType>
                SubtractInternalImagesType;
        typename SubtractInternalImagesType::Pointer vBGResidualImageFilter = SubtractInternalImagesType::New();
        vBGResidualImageFilter->SetInput1(this->m_DataImage);
        vBGResidualImageFilter->SetInput2(vModelImageSubtractionFilter->GetOutput());
        vBGResidualImageFilter->Update();

        // Time vs. Memory tradeoff:
        // Memory efficient: iterate the label image: for all new seed points (new label found),
        // iterate through the label usinga floodfill iterator. While iterating,
        // read out the data and model image. Put the quotient of those into an
        // 'new' array. Sort the array, read out the median and delete the array.
        //
        // Time efficient: iterate the label image. Store all quotients found in
        // an array corresponding to the label. This is another 32-bit copy of
        // the image.


        /// Set up a map datastructure that maps from labels to arrays of data
        /// values. These arrays will be sorted to read out the median.
        std::map <LabelAbsPixelType, unsigned int> vLabelCounter;
        std::map <LabelAbsPixelType, float> vIntensitySum;

        using namespace boost::accumulators;
        typedef accumulator_set<InternalPixelType, stats<tag::median> > AccType;         
        typedef std::map <LabelAbsPixelType, AccType> LabelToAccMapType;
        LabelToAccMapType vScalings3;
        
        /// For all the active labels, create an entry in the map and initialize
        /// an array as the corresponding value.
        CountStatisticsIteratorType vActiveLabelsIt =
                this->m_Count.begin();
        CountStatisticsIteratorType vActiveLabelsItEnd =
                this->m_Count.end();
        for (; vActiveLabelsIt != vActiveLabelsItEnd; ++vActiveLabelsIt) {
            LabelAbsPixelType vLabelAbs = vActiveLabelsIt->first;
//            if (vLabelAbs == static_cast<unsigned int> (m_ForbiddenRegionLabel)) {
//                continue;
//            }
            vScalings3[vLabelAbs] = AccType();
            vLabelCounter[vLabelAbs] = 0;
            vIntensitySum[vLabelAbs] = 0.0f;
        }

        ImageRegionConstIterator<InternalImageType> vScaleIt(vDivisionFilter->GetOutput(),
                vDivisionFilter->GetOutput()->GetBufferedRegion());
        ImageRegionConstIterator<LabelImageType> vLabelIt(this->m_LabelImage,
                this->m_LabelImage->GetBufferedRegion());
        ImageRegionConstIterator<InternalImageType> vBGIt(vBGResidualImageFilter->GetOutput(),
                vBGResidualImageFilter->GetOutput()->GetBufferedRegion());


        for (vScaleIt.GoToBegin(), vLabelIt.GoToBegin(), vBGIt.GoToBegin();
                !vLabelIt.IsAtEnd();
                ++vScaleIt, ++vLabelIt, ++vBGIt) {
            LabelAbsPixelType vLabelAbs = abs(vLabelIt.Get());
//            if (vLabelAbs == static_cast<unsigned int> (m_ForbiddenRegionLabel)) {
//                continue;
//            }
            unsigned int vC = vLabelCounter[vLabelAbs]++;

            if (vLabelAbs == 0) {
                vScalings3[vLabelAbs](vBGIt.Get());
            } else {
                vScalings3[vLabelAbs](vScaleIt.Get());
            }
        }

        /// For all the labels (except the BG ?) sort the scalar factors for all
        /// the pixel. The median is in the middle of the sorted list.
        /// TODO: check if fitting is necessary for the BG.
        /// TODO: Depending on the size of the region, sorting takes too long and
        ///       a median of medians algorithm (O(N)) could provide a good
        ///       approximation of the median.
        typename LabelToAccMapType::iterator vLSCIt = vScalings3.begin();
        typename LabelToAccMapType::iterator vLSCItEnd = vScalings3.end();

        for (; vLSCIt != vLSCItEnd; ++vLSCIt) {
            LabelAbsPixelType vLabelAbs = vLSCIt->first;
            float vMedian = 0;
            if (this->m_Count[vLabelAbs] > 2) {
                vMedian = median(vScalings3[vLabelAbs]); // median from boost
            } else {
                if(this->m_Count[vLabelAbs] > 0){
                    vMedian = this->m_Sums[vLabelAbs]/this->m_Count[vLabelAbs];
                }
            }

            /// Correct the old intensity values.
            if (vLabelAbs == 0) {
                if (vMedian < 0) { // TODO: only accept unsigned data types
                    vMedian = 0;
                }
                m_Intensities[vLabelAbs] = vMedian;
            } else {
                m_Intensities[vLabelAbs] =
                        (m_Intensities[vLabelAbs] - vOldBG) * vMedian + m_Intensities[0];
            }
        }

        // The model image has to be renewed as well to match the new statistic values:
        void * vimg = &m_DeconvolutionModelImage;
        GenerateModelImage(vimg, this->m_LabelImage, this->m_DataImage, m_PSF, &m_Intensities);
        
    }

    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void 
    RCPiecewiseConstantDeconvolutionEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
    ::GenerateModelImage(void* aPointerToResultImage, 
            LabelImageConstPointerType aLabelImage, 
            InternalImageConstPointerType aDataImage, 
            InternalImagePointerType aPSFImage, 
            SumStatisticsMapType* aIntensityValues) {
        
        typename InternalImageType::Pointer vModelImage = InternalImageType::New();

        vModelImage->SetRequestedRegion(aDataImage->GetRequestedRegion());
        vModelImage->SetBufferedRegion(aDataImage->GetBufferedRegion());
        vModelImage->SetLargestPossibleRegion(aDataImage->GetLargestPossibleRegion());
        vModelImage->Allocate();
        vModelImage->SetSpacing(this->m_DataImage->GetSpacing());
        ImageRegionConstIterator<LabelImageType> vLabelIt(aLabelImage,
                aLabelImage->GetBufferedRegion());
        ImageRegionIterator<InternalImageType> vConvIt(vModelImage,
                aLabelImage->GetBufferedRegion());
        
        for (vLabelIt.GoToBegin(), vConvIt.GoToBegin();
                !vLabelIt.IsAtEnd();
                ++vLabelIt, ++vConvIt) {
            LabelPixelType vLabel = vLabelIt.Get();
//            if (vLabel == m_ForbiddenRegionLabel) {
//                vLabel = 0; // Set Background value ??
//            }
            vConvIt.Set((*aIntensityValues)[abs(vLabel)]);
            
        }
        
        typedef FFTConvolutionImageFilter<InternalImageType, InternalImageType> ConvolutionFilterType;
        typename ConvolutionFilterType::Pointer vConvolutionFilter = ConvolutionFilterType::New();
        vConvolutionFilter->SetKernelImage(aPSFImage);
        vConvolutionFilter->SetInput(vModelImage);
        vConvolutionFilter->Update();
        
        
        typename InternalImageType::Pointer* vReturnPointer =
        (typename InternalImageType::Pointer*) aPointerToResultImage;
        *vReturnPointer = vConvolutionFilter->GetOutput();
    }    
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void 
    RCPiecewiseConstantDeconvolutionEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
    ::UpdateModelImage(IndexType aIndex,
    LabelAbsPixelType aFromLabel,
    LabelAbsPixelType aToLabel) {

        if(m_OptimizationMode) {
            return; // the model image will be regenerated from scratch, speedup.
        }
        
        typedef typename LabelImageType::RegionType RegionType;
        RegionType vRegion;
        typename LabelImageType::SizeType vRegionSize = m_PSF->GetLargestPossibleRegion().GetSize();

        typename LabelImageType::OffsetType vOffset;
        vOffset.Fill(vRegionSize[0] / 2); // we assume the region to be a hypercube.
        vRegion.SetSize(vRegionSize);
        vRegion.SetIndex(aIndex - vOffset);

        /// Move the PSF window such that the center is on aIndex
        RegionType vPSFRegion = m_PSF->GetLargestPossibleRegion().GetSize();
        vPSFRegion.SetIndex(vRegion.GetIndex());
        m_PSF->SetBufferedRegion(vPSFRegion);

        /// After the cropping at the data-image boundaries, vRegion will be the
        /// region to treat in the data-image space.
        vRegion.Crop(m_DeconvolutionModelImage->GetBufferedRegion());

        /// Iterate through the region and subtract the PSF from the conv image.
        ImageRegionIterator<InternalImageType> vIt(m_DeconvolutionModelImage, vRegion);
        ImageRegionConstIterator<InternalImageType> vPSFIt(m_PSF, vRegion);

        // If the new label is the BG then subtract a PSf from the model image.
        // Else remove the PSF from the old and add a PSF from the new region.
        InternalPixelType vIntensity_From = m_Intensities[aFromLabel];
        InternalPixelType vIntensity_To = m_Intensities[aToLabel];
        for (vPSFIt.GoToBegin(), vIt.GoToBegin();
                !vPSFIt.IsAtEnd();
                ++vPSFIt, ++vIt) {
            vIt.Set(vIt.Get() + (vIntensity_To - vIntensity_From) * vPSFIt.Get());
        }
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>    
    typename RCPiecewiseConstantDeconvolutionEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>::
    EnergyDifferenceType
    RCPiecewiseConstantDeconvolutionEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
    ::CalculateTotalEnergy(){
        EnergyDifferenceType vEnergy = 0;
        /* Outline:
         * Generate the 'ideal image' J:
         * - Set a BG valueCalcFullEnerg
         * - Iterate the label image and put the corresponding intensity value
         *   to the corresponding iterator position.
         * - Convolve with the point-spread-function
         */
        typename InternalImageType::Pointer vIdealImage;
        void * vimg = &vIdealImage;
        GenerateModelImage(vimg, this->m_LabelImage, this->m_DataImage, this->m_PSF, &m_Intensities);

        ImageRegionConstIteratorWithIndex<DataImageType> vDataIt(this->m_DataImage,
                this->m_DataImage->GetBufferedRegion());
        ImageRegionConstIterator<InternalImageType> vIdealIt(vIdealImage,
                vIdealImage->GetBufferedRegion());
        ImageRegionConstIterator<LabelImageType> vLabelIt(this->m_LabelImage,
                this->m_LabelImage->GetBufferedRegion());

        for (vIdealIt.GoToBegin(), vDataIt.GoToBegin(), vLabelIt.GoToBegin();
                !vIdealIt.IsAtEnd();
                ++vIdealIt, ++vDataIt, ++vLabelIt) {
            if (this->m_Coefficient && vLabelIt.Get() != NumericTraits<LabelPixelType>::max()) {
                vEnergy += this->m_Coefficient * (vIdealIt.Get() - vDataIt.Get()) *
                        (vIdealIt.Get() - vDataIt.Get());
            }
        }
        return vEnergy;
    }

    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>    
    void
    RCPiecewiseConstantDeconvolutionEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
    ::GenerateReconstructedImage(void* aPointerToResultImage){
        
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
            InternalPixelType vAbsLabel = abs(vPixelIt.Get());
            IndexType vCenterIndex = vPixelIt.GetIndex();
            vReconstImage->SetPixel(vCenterIndex, m_Intensities[vAbsLabel]);
        }
        
        /* set the smart pointer */
        typename ReconstructedImageType::Pointer* vReturnPointer =
        (typename ReconstructedImageType::Pointer*) aPointerToResultImage;
        *vReturnPointer = vReconstImage;
        return;
    }

} // end namespace itk
