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
#ifndef __itkSphereBitmapImageSource_hxx
#define __itkSphereBitmapImageSource_hxx

#include "itkSphereBitmapImageSource.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkProgressReporter.h"
#include "itkObjectFactory.h"

// template from:
//#include "itkGaussianImageSource.h"

namespace itk {

    template <class TOutputImage>
    SphereBitmapImageSource<TOutputImage>
    ::SphereBitmapImageSource() {
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

        m_Radius = 10;

        m_BackgroundValue = static_cast<OutputImagePixelType> (0);
        m_ForegroundValue = static_cast<OutputImagePixelType> (1);

    }

    template <class TOutputImage>
    SphereBitmapImageSource<TOutputImage>
    ::~SphereBitmapImageSource() {
    }

    template <class TOutputImage>
    void
    SphereBitmapImageSource<TOutputImage>
    ::PrintSelf(std::ostream& os, Indent indent) const {
        Superclass::PrintSelf(os, indent);

        os << indent << "ForegroundValue: " << m_ForegroundValue << std::endl;
        os << indent << "BackgroundValue: " << m_BackgroundValue << std::endl;

    }

    //----------------------------------------------------------------------------

    template <typename TOutputImage>
    void
    SphereBitmapImageSource<TOutputImage>
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
    SphereBitmapImageSource<TOutputImage>
    ::GenerateData() {
        typename TOutputImage::Pointer outputPtr = this->GetOutput();

        // allocate the output buffer
        outputPtr->SetBufferedRegion(outputPtr->GetRequestedRegion());
        outputPtr->Allocate();

        // Create an iterator that will walk the output region
        typedef ImageRegionIteratorWithIndex<TOutputImage> OutputIterator;
        OutputIterator outIt = OutputIterator(outputPtr,
                outputPtr->GetRequestedRegion());

        //        ProgressReporter progress(this, 0,
        //                outputPtr->GetRequestedRegion()
        //                .GetNumberOfPixels());
        unsigned int vDim = TOutputImage::ImageDimension;

        // TODO: see itkFlatStructuringElement: use an 
        // EllipsoidInteriorExteriorSpatialFunction with a flood-fill iterator. 
        // OR: dilate a point with FlatStructuringElement::Ball.
        for (outIt.GoToBegin(); !outIt.IsAtEnd(); ++outIt) {
            IndexType vIndex = outIt.GetIndex();

            float vHypEllipse = 0;
            for (unsigned int vD = 0; vD < vDim; vD++) {
                vHypEllipse +=
                        (static_cast<double> (vIndex[vD]) - static_cast<double> (m_Size[vD]-1) / 2.0) *
                        (static_cast<double> (vIndex[vD]) - static_cast<double> (m_Size[vD]-1) / 2.0) /
                        (m_Radius[vD] * m_Radius[vD]);
            }

            if (vHypEllipse <= 1.0f) {
                outIt.Set(static_cast<OutputImagePixelType> (m_ForegroundValue));
            } else {
                outIt.Set(static_cast<OutputImagePixelType> (m_BackgroundValue));
            }
        }
        
    }

    template<typename TOutputImage>
    void
    SphereBitmapImageSource<TOutputImage>
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
    SphereBitmapImageSource<TOutputImage>
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
    SphereBitmapImageSource<TOutputImage>
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
    SphereBitmapImageSource<TOutputImage>
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
    SphereBitmapImageSource<TOutputImage>
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
    SphereBitmapImageSource<TOutputImage>
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
