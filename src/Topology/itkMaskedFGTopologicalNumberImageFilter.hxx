/*=========================================================================
 *
 * This code has been downloaded from 
 * http://www.insight-journal.org/browse/publication/120.
 * 
 * The code was written and contributed by Lamy Julien, 2006.
 * 
 *=========================================================================*/
#ifndef __ITKMASKEDFGTOPOLOGICALNUMBERIMAGEFILTER_HXX_
#define __ITKMASKEDFGTOPOLOGICALNUMBERIMAGEFILTER_HXX_

#include "itkMaskedFGTopologicalNumberImageFilter.h"
#include "itkTopologicalNumberImageFunction.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"

namespace itk
{

    template <class TInputImage, class TOutputImage>
            MaskedFGTopologicalNumberImageFilter<TInputImage, TOutputImage>
            ::MaskedFGTopologicalNumberImageFilter() {

        m_TopologicalNumberFunction = TopologicalNumberImageFunction<
                InputImageType, ForegroundConnectivityType>::New();

//        m_BackgroundValue = static_cast<OutputPixelType>(0);
//        m_ForegroundValue = static_cast<OutputPixelType>(1);
    }

    template <class TInputImage, class TOutputImage>
            MaskedFGTopologicalNumberImageFilter<TInputImage, TOutputImage>
            ::~MaskedFGTopologicalNumberImageFilter() {
    }

    //TODO:
//    template <class TInputImage, class TOutputImage>
//    void MaskedFGTopologicalNumberImageFilter<TInputImage, TOutputImage>
//    ::ThreadedGenerateData(const OutputImageRegionType&
//    outputRegionForThread, int threadId) {
//
//    }
//

    template <class TInputImage, class TOutputImage>
    void MaskedFGTopologicalNumberImageFilter<TInputImage, TOutputImage>
    ::GenerateData() {

        /**
         * Set up the regions and allocate the output-image
         */
        this->GetOutput()->SetRequestedRegion(this->GetInput()->GetRequestedRegion());
        this->GetOutput()->SetBufferedRegion(this->GetInput()->GetBufferedRegion());
        this->GetOutput()->SetLargestPossibleRegion(this->GetInput()->GetLargestPossibleRegion());
        this->GetOutput()->Allocate();

        typedef itk::ImageRegionConstIteratorWithIndex<TInputImage> InputImageIteratorType;
        InputImageIteratorType vInputIt(this->GetInput(),
                this->GetInput()->GetBufferedRegion());

        typedef itk::ImageRegionIterator<TOutputImage> OutputImageIteratorType;
        OutputImageIteratorType vOutputIt(this->GetOutput(),
                this->GetOutput()->GetBufferedRegion());


        m_TopologicalNumberFunction->SetInputImage(this->GetInput());

        for(vInputIt.GoToBegin(), vOutputIt.GoToBegin(); !vInputIt.IsAtEnd();
                ++vInputIt, ++vOutputIt) {

            if (vInputIt.Get() != NumericTraits<InputPixelType>::Zero) {
                std::pair<unsigned int, unsigned int> vTopoNbPair =
                        m_TopologicalNumberFunction->EvaluateAtIndex(vInputIt.GetIndex());
                if(vTopoNbPair.first == 2 && vTopoNbPair.second == 1) {
                    vOutputIt.Set(2); //curve point
                } else if(vTopoNbPair.first == 1 && vTopoNbPair.second == 2) {
                    vOutputIt.Set(3); // surface point
//                    std::cout << "Warning: surface point detected at " << vInputIt.GetIndex() << std::endl;
                } else if(vTopoNbPair.first > 1 && vTopoNbPair.second == 1) {
                    vOutputIt.Set(4); // junction point
                } else if(vTopoNbPair.first == 1 && vTopoNbPair.second > 2) {
                    vOutputIt.Set(5);
                } else if (vTopoNbPair.first > 2 && vTopoNbPair.second == 1) {
                    vOutputIt.Set(6);
                } else if(vTopoNbPair.first == 1 && vTopoNbPair.second == 1) {
                    std::cout << "Simple point detected at " << vInputIt.GetIndex() << std::endl;
                    vOutputIt.Set(1);
                }
            }
        }
    }

    template <class TInputImage, class TOutputImage>
            void MaskedFGTopologicalNumberImageFilter<TInputImage, TOutputImage>
            ::PrintSelf(std::ostream& aOS, Indent aIndent) const {
//        aOS << aIndent << "ForegroundValue: " << m_ForegroundValue;
//        aOS << aIndent << "BackgroundValue: " << m_BackgroundValue;

    }

}//end namespace itk

#endif //__ITKMASKEDFGTOPOLOGICALNUMBERIMAGEFILTER_HXX_
