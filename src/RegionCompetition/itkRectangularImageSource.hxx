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
#ifndef __itkRectangularImageSource_hxx
#define __itkRectangularImageSource_hxx

#include "itkRectangularImageSource.h"
#include "itkImageRegionIterator.h"
#include "itkProgressReporter.h"
#include "itkObjectFactory.h"

// template from:
//#include "itkGaussianImageSource.h"

namespace itk {

    template <class TOutputImage>
    RectangularImageSource<TOutputImage>
    ::RectangularImageSource() {
        SizeType vSize;
        IndexType vIndex;

        //Initial image is 64 wide in each direction.
        for (unsigned int i = 0; i < TOutputImage::GetImageDimension(); i++) {
            m_Size[i] = 64;
            m_Spacing[i] = 1.0;
            m_Origin[i] = 0.0;
            vSize[i] = 0;
            vIndex[i] = 0;
        }
        m_Direction.SetIdentity();

        m_ForegroundRegion.SetIndex(vIndex);
        m_ForegroundRegion.SetSize(vSize);

        m_BackgroundValue = static_cast<OutputImagePixelType> (0);
        m_ForegroundValue = static_cast<OutputImagePixelType> (1);

    }

    template <class TOutputImage>
    RectangularImageSource<TOutputImage>
    ::~RectangularImageSource() {
    }

    template <class TOutputImage>
    void
    RectangularImageSource<TOutputImage>
    ::PrintSelf(std::ostream& os, Indent indent) const {
        Superclass::PrintSelf(os, indent);

        os << indent << "ForegroundValue: " << m_ForegroundValue << std::endl;
        os << indent << "BackgroundValue: " << m_BackgroundValue << std::endl;
        os << indent << "ForegroundRegion: " << m_ForegroundRegion << std::endl;

    }

    //----------------------------------------------------------------------------

    template <typename TOutputImage>
    void
    RectangularImageSource<TOutputImage>
    ::GenerateOutputInformation() {
        TOutputImage *output;
        typename TOutputImage::IndexType index = {
            {0}
        };
        typename TOutputImage::SizeType size = {
            {0}
        };
        size.SetSize(m_Size);

        output = this->GetOutput(0);

        typename TOutputImage::RegionType largestPossibleRegion;
        largestPossibleRegion.SetSize(size);
        largestPossibleRegion.SetIndex(index);
        output->SetLargestPossibleRegion(largestPossibleRegion);

        output->SetSpacing(m_Spacing);
        output->SetOrigin(m_Origin);
        output->SetDirection(m_Direction);
    }

    template <typename TOutputImage>
    void
    RectangularImageSource<TOutputImage>
    ::GenerateData() {
        typename TOutputImage::Pointer outputPtr = this->GetOutput();

        // allocate the output buffer
        outputPtr->SetBufferedRegion(outputPtr->GetRequestedRegion());
        outputPtr->Allocate();

        // Create an iterator that will walk the output region
        typedef ImageRegionIterator<TOutputImage> OutputIterator;
        OutputIterator outIt = OutputIterator(outputPtr,
                outputPtr->GetRequestedRegion());

        //        ProgressReporter progress(this, 0,
        //                outputPtr->GetRequestedRegion()
        //                .GetNumberOfPixels());

        // Walk the output image, evaluating the spatial function at each pixel
        for (outIt.GoToBegin(); !outIt.IsAtEnd(); ++outIt) {
            // Set the pixel value to the background value
            outIt.Set(static_cast<OutputImagePixelType> (m_BackgroundValue));
        }

        OutputIterator vRectIt = OutputIterator(outputPtr, m_ForegroundRegion);
        for (vRectIt.GoToBegin(); !vRectIt.IsAtEnd(); ++vRectIt) {
            // Set the pixel inside the foreground region to the foreground value
            vRectIt.Set(static_cast<OutputImagePixelType> (m_ForegroundValue));
        }

    }

    template<typename TOutputImage>
    void
    RectangularImageSource<TOutputImage>
    ::SetSpacing(const float* spacing) {
        unsigned int i;
        for (i = 0; i < TOutputImage::ImageDimension; i++) {
            if ((double) spacing[i] != m_Spacing[i]) {
                break;
            }
        }
        if (i < TOutputImage::ImageDimension) {
            for (i = 0; i < TOutputImage::ImageDimension; i++) {
                m_Spacing[i] = spacing[i];
            }
            this->Modified();
        }
    }

    template<typename TOutputImage>
    void
    RectangularImageSource<TOutputImage>
    ::SetSpacing(const double* spacing) {
        unsigned int i;
        for (i = 0; i < TOutputImage::ImageDimension; i++) {
            if (spacing[i] != m_Spacing[i]) {
                break;
            }
        }
        if (i < TOutputImage::ImageDimension) {
            for (i = 0; i < TOutputImage::ImageDimension; i++) {
                m_Spacing[i] = spacing[i];
            }
            this->Modified();
        }
    }

    template<typename TOutputImage>
    void
    RectangularImageSource<TOutputImage>
    ::SetOrigin(const float* origin) {
        unsigned int i;
        for (i = 0; i < TOutputImage::ImageDimension; i++) {
            if ((double) origin[i] != m_Origin[i]) {
                break;
            }
        }
        if (i < TOutputImage::ImageDimension) {
            for (i = 0; i < TOutputImage::ImageDimension; i++) {
                m_Origin[i] = origin[i];
            }
            this->Modified();
        }
    }

    template<typename TOutputImage>
    void
    RectangularImageSource<TOutputImage>
    ::SetOrigin(const double* origin) {
        unsigned int i;
        for (i = 0; i < TOutputImage::ImageDimension; i++) {
            if (origin[i] != m_Origin[i]) {
                break;
            }
        }
        if (i < TOutputImage::ImageDimension) {
            for (i = 0; i < TOutputImage::ImageDimension; i++) {
                m_Origin[i] = origin[i];
            }
            this->Modified();
        }
    }

    template<typename TOutputImage>
    void
    RectangularImageSource<TOutputImage>
    ::SetSize(const SizeValueType * size) {
        unsigned int i;
        for (i = 0; i < TOutputImage::ImageDimension; i++) {
            if (size[i] != m_Size[i]) {
                break;
            }
        }
        if (i < TOutputImage::ImageDimension) {
            for (i = 0; i < TOutputImage::ImageDimension; i++) {
                m_Size[i] = size[i];
            }
            this->Modified();
        }
    }

    template<typename TOutputImage>
    void
    RectangularImageSource<TOutputImage>
    ::SetSize(const SizeType size) {
        unsigned int i;
        for (i = 0; i < TOutputImage::ImageDimension; i++) {
            if (size[i] != m_Size[i]) {
                break;
            }
        }
        if (i < TOutputImage::ImageDimension) {
            for (i = 0; i < TOutputImage::ImageDimension; i++) {
                m_Size[i] = size[i];
            }
            this->Modified();
        }
    }


} // end namespace itk

#endif
