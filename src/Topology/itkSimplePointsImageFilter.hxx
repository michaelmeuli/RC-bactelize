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

#ifndef __ITKSIMPLEPOINTSIMAGEFILTER_CXX_
#define __ITKSIMPLEPOINTSIMAGEFILTER_CXX_

#include "itkMaskedFGTopologicalNumberImageFilter.h"


#include "itkSimplePointsImageFilter.h"
#include "itkSimplicityByTopologicalNumbersImageFunction.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
namespace itk
{

    template <class TInputImage, class TOutputImage>
            SimplePointsImageFilter<TInputImage, TOutputImage>
            ::SimplePointsImageFilter() {
        m_SimplicityCriterion = SimplicityByTopologicalNumbersImageFunction<
                InputImageType, ForegroundConnectivityType>::New();

        m_BackgroundValue = static_cast<OutputPixelType>(0);
        m_ForegroundValue = static_cast<OutputPixelType>(1);
    }

    template <class TInputImage, class TOutputImage>
            SimplePointsImageFilter<TInputImage, TOutputImage>
            ::~SimplePointsImageFilter() {
    }

//    template <class TInputImage, class TOutputImage>
//    void SimplePointsImageFilter<TInputImage, TOutputImage>
//    ::ThreadedGenerateData(const OutputImageRegionType&
//    outputRegionForThread, int threadId) {
//
//        typedef itk::ImageRegionConstIteratorWithIndex<TInputImage> InputImageIteratorType;
//        InputImageIteratorType vInputIt(this->GetInput(),
//                outputRegionForThread);
//
//        typedef itk::ImageRegionIterator<TOutputImage> OutputImageIteratorType;
//        OutputImageIteratorType vOutputIt(this->GetOutput(),
//                outputRegionForThread);
//
//
//        m_SimplicityCriterion->SetInputImage(this->GetInput());
//
//        for (vInputIt.GoToBegin(), vOutputIt.GoToBegin(); !vInputIt.IsAtEnd();
//                ++vInputIt, ++vOutputIt) {
//            InputImageIndexType vIndex = vInputIt.GetIndex();
//
//            if (m_SimplicityCriterion->EvaluateAtIndex(vIndex)) {
//                vOutputIt.Set(m_ForegroundValue);
//            } else {
//                vOutputIt.Set(m_BackgroundValue);
//            }
//
//        }
//    }
//
    template <class TInputImage, class TOutputImage>
            void SimplePointsImageFilter<TInputImage, TOutputImage>
            ::GenerateData() {

        /**
         * Set up the regions and allocate the output-image
         */
        this->GetOutput()->SetRequestedRegion(this->GetInput()->GetRequestedRegion());
        this->GetOutput()->SetBufferedRegion(this->GetInput()->GetBufferedRegion());
        this->GetOutput()->SetLargestPossibleRegion(this->GetInput()->GetLargestPossibleRegion());
        this->GetOutput()->Allocate();

        typename InputImageType::RegionType vRegionWithoutBoundary = this->GetInput()->GetBufferedRegion();
        typename InputImageType::SizeType vReducedSize = vRegionWithoutBoundary.GetSize();
        typename InputImageType::IndexType vOffsetIndex = vRegionWithoutBoundary.GetIndex();
        for(unsigned int vD = 0; vD < TInputImage::ImageDimension; vD++) {
            //TODO unsafe! region could of negative size and outside of the buffered region
            vReducedSize[vD] = vReducedSize[vD] - 2;
            vOffsetIndex[vD] = vOffsetIndex[vD] + 1;
        }
        vRegionWithoutBoundary.SetSize(vReducedSize);
        vRegionWithoutBoundary.SetIndex(vOffsetIndex);
        
        typedef itk::ImageRegionConstIteratorWithIndex<TInputImage> InputImageIteratorType;
        InputImageIteratorType vInputIt(this->GetInput(),
                vRegionWithoutBoundary);
//                this->GetInput()->GetBufferedRegion());

        typedef itk::ImageRegionIterator<TOutputImage> OutputImageIteratorType;
        OutputImageIteratorType vOutputIt(this->GetOutput(),
                vRegionWithoutBoundary);
//                this->GetOutput()->GetBufferedRegion());


        m_SimplicityCriterion->SetInputImage(this->GetInput());

        for(vInputIt.GoToBegin(), vOutputIt.GoToBegin(); !vInputIt.IsAtEnd();
                ++vInputIt, ++vOutputIt) {
            InputImageIndexType vIndex = vInputIt.GetIndex();

            if (m_SimplicityCriterion->EvaluateAtIndex(vIndex)) {
                vOutputIt.Set(m_ForegroundValue);
            } else {
                vOutputIt.Set(m_BackgroundValue);
            }

        }
    }

    template <class TInputImage, class TOutputImage>
            void SimplePointsImageFilter<TInputImage, TOutputImage>
            ::PrintSelf(std::ostream& aOS, Indent aIndent) const {
        aOS << aIndent << "ForegroundValue: " << m_ForegroundValue;
        aOS << aIndent << "BackgroundValue: " << m_BackgroundValue;

    }

}//end namespace itk

#endif
