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
#ifndef __itkRectangularImageSource_h
#define __itkRectangularImageSource_h

#include "itkImageSource.h"
#include "itkFixedArray.h"
#include "itkSize.h"

namespace itk {

    /** \class RectangularImageSource
     * \brief Generate an n-dimensional image of a Foreground.
     *
     *
     * The output image may be of any dimension.
     *
     * \ingroup DataSources
     */
    template <typename TOutputImage>
    class ITK_EXPORT RectangularImageSource : public ImageSource<TOutputImage> {
    public:
        /** Standard class typedefs. */
        typedef RectangularImageSource Self;
        typedef ImageSource<TOutputImage> Superclass;
        typedef SmartPointer<Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;

        /** Typedef for the output image PixelType. */
        typedef typename TOutputImage::PixelType OutputImagePixelType;

        /** Typedef to describe the output image region type. */
        typedef typename TOutputImage::RegionType OutputImageRegionType;

        /** Spacing typedef support.  Spacing holds the size of a pixel.  The
         * spacing is the geometric distance between image samples. */
        typedef typename TOutputImage::SpacingType SpacingType;

        /** Origin typedef support.  The origin is the geometric coordinates
         * of the index (0,0). */
        typedef typename TOutputImage::PointType PointType;

        /** Index typedef support */
        typedef typename TOutputImage::IndexType IndexType;

        /** Direction typedef support.  The direction is the direction
         * cosines of the image. */
        typedef typename TOutputImage::DirectionType DirectionType;

        /** Dimensionality of the output image */
        itkStaticConstMacro(NDimensions, unsigned int, TOutputImage::ImageDimension);

        /** Type used to store gaussian parameters. */
        typedef FixedArray<double, itkGetStaticConstMacro(NDimensions) > ArrayType;

        /** Size type matches that used for images */
        typedef typename TOutputImage::SizeType SizeType;
        typedef typename TOutputImage::SizeValueType SizeValueType;

        /** Run-time type information (and related methods). */
        itkTypeMacro(RectangularImageSource, ImageSource);

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** Specify the size of the output image. */
        virtual void SetSize(const SizeValueType * values);

        /** Specify the size of the output image. */
        virtual void SetSize(const SizeType values);

        /** Get the size of the output image. */
        itkGetVectorMacro(Size, const SizeValueType, NDimensions);

        /** Specify the spacing of the output image. */
        itkSetMacro(Spacing, SpacingType);
        virtual void SetSpacing(const float* values);
        virtual void SetSpacing(const double* values);

        /** Get the spacing of the output image. */
        itkGetConstReferenceMacro(Spacing, SpacingType);

        /** Specify the origin of the output image. */
        itkSetMacro(Origin, PointType);
        virtual void SetOrigin(const float* values);
        virtual void SetOrigin(const double* values);

        /** Get the origin of the output image. */
        itkGetConstReferenceMacro(Origin, PointType);

        /** Specify the direction of the output image. */
        itkSetMacro(Direction, DirectionType);
        itkGetConstReferenceMacro(Direction, DirectionType);

        /** Gets and sets for gaussian parameters */
        itkSetMacro(ForegroundRegion, OutputImageRegionType);
        itkGetConstReferenceMacro(ForegroundRegion, OutputImageRegionType);
        itkSetMacro(ForegroundValue, OutputImagePixelType);
        itkGetConstReferenceMacro(ForegroundValue, OutputImagePixelType);
        itkSetMacro(BackgroundValue, OutputImagePixelType);
        itkGetConstReferenceMacro(BackgroundValue, OutputImagePixelType);


    protected:
        RectangularImageSource();
        ~RectangularImageSource();
        void PrintSelf(std::ostream& os, Indent indent) const;
        void GenerateData();
        virtual void GenerateOutputInformation();

    private:
        RectangularImageSource(const RectangularImageSource&); //purposely not implemented
        void operator=(const RectangularImageSource&); //purposely not implemented

        SizeValueType m_Size[NDimensions]; //size of the output image
        SpacingType m_Spacing; //spacing
        PointType m_Origin; //origin
        DirectionType m_Direction; // direction

        /** Parameters for the Rectangle. */
        OutputImageRegionType m_ForegroundRegion;
        OutputImagePixelType m_ForegroundValue;
        OutputImagePixelType m_BackgroundValue;

    };

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRectangularImageSource.hxx"
#endif

#endif
