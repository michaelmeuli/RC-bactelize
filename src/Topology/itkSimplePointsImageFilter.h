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

#ifndef _ITKSIMPLEPOINTSIMAGEFILTER_H
#define	_ITKSIMPLEPOINTSIMAGEFILTER_H

#include "itkImageToImageFilter.h"
#include "itkImageFunction.h"
#include "itkImage.h"
#include "itkConnectivity.h"

/// TODO: template with foreground connectivity type (default 0)
namespace itk {

    template <class TInputImage, class TOutputImage >
    class ITK_EXPORT SimplePointsImageFilter :
    public ImageToImageFilter<TInputImage, TOutputImage> {

    public:
        /** Standard class typedefs         */
        typedef SimplePointsImageFilter Self;
        typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
        typedef SmartPointer<Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** Image Types */
        typedef TInputImage InputImageType;
        typedef TOutputImage OutputImageType;

        /** Typedef to describe the output and input image region types. */
        typedef typename InputImageType::RegionType InputImageRegionType;
        typedef typename OutputImageType::RegionType OutputImageRegionType;

        /** Typedef to describe the type of pixel. */
        typedef typename InputImageType::PixelType InputPixelType;
        typedef typename OutputImageType::PixelType OutputPixelType;

        /** Typedef to describe the output and input image index and size types. */
        typedef typename InputImageType::IndexType InputImageIndexType;
        typedef typename InputImageType::SizeType InputImageSizeType;
        typedef typename InputImageType::OffsetType InputImageOffsetType;
        typedef typename OutputImageType::IndexType OutputImageIndexType;
        typedef typename OutputImageType::SizeType OutputImageSizeType;
        typedef typename OutputImageType::OffsetType OutputImageOffsetType;

    protected:
        SimplePointsImageFilter();
        virtual ~SimplePointsImageFilter();
        void PrintSelf(std::ostream& aOS, Indent aIndent) const;

// TODO: doesnt work, most probably because of errors at the boundary of each region:
//        void ThreadedGenerateData(const OutputImageRegionType&
//        		aOutputRegionForThread, int aThreadId);

        void GenerateData();

    private:
        typedef typename itk::Connectivity<InputImageType::ImageDimension, InputImageType::ImageDimension-1> ForegroundConnectivityType;
        typename ImageFunction<InputImageType, bool >::Pointer m_SimplicityCriterion;
        OutputPixelType m_ForegroundValue;
        OutputPixelType m_BackgroundValue;
    public:
        itkSetMacro(ForegroundValue, OutputPixelType);
        itkSetMacro(BackgroundValue, OutputPixelType);
        itkGetConstMacro(ForegroundValue, OutputPixelType);
        itkGetConstMacro(BackgroundValue, OutputPixelType);

    };

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSimplePointsImageFilter.hxx"
#endif

#endif /* _ITKSIMPLEPOINTSIMAGEFILTER_H */
