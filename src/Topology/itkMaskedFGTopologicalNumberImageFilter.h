/*=========================================================================
 *
 * This code has been downloaded from 
 * http://www.insight-journal.org/browse/publication/120.
 * 
 * The code was originally written and contributed by Lamy Julien, 2006.
 *=========================================================================*/
#ifndef _ITKMASKEDFGTOPOLOGICALNUMBERIMAGEFILTER_H
#define	_ITKMASKEDFGTOPOLOGICALNUMBERIMAGEFILTER_H

#include "itkImageToImageFilter.h"
#include "itkImageFunction.h"
#include "itkImage.h"
#include "itkConnectivity.h"

/// TODO: template with foreground connectivity type (default DIM-1)
namespace itk {

    template <class TInputImage, class TOutputImage >
    class ITK_EXPORT MaskedFGTopologicalNumberImageFilter :
    public ImageToImageFilter<TInputImage, TOutputImage> {

    public:
        /** Standard class typedefs         */
        typedef MaskedFGTopologicalNumberImageFilter Self;
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
        MaskedFGTopologicalNumberImageFilter();
        virtual ~MaskedFGTopologicalNumberImageFilter();
        void PrintSelf(std::ostream& aOS, Indent aIndent) const;

// TODO: doesnt work, most probably because of errors at the boundary of each region:
//        void ThreadedGenerateData(const OutputImageRegionType&
//        		aOutputRegionForThread, int aThreadId);

        void GenerateData();

    private:
        typedef typename itk::Connectivity<InputImageType::ImageDimension, InputImageType::ImageDimension-1> ForegroundConnectivityType;
        typename ImageFunction<InputImageType, std::pair<unsigned int, unsigned int> >::Pointer m_TopologicalNumberFunction;
//        OutputPixelType m_ForegroundValue;
//        OutputPixelType m_BackgroundValue;
    public:
//        itkSetMacro(ForegroundValue, OutputPixelType);
//        itkSetMacro(BackgroundValue, OutputPixelType);
//        itkGetConstMacro(ForegroundValue, OutputPixelType);
//        itkGetConstMacro(BackgroundValue, OutputPixelType);

    };

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMaskedFGTopologicalNumberImageFilter.hxx"
#endif

#endif /* _ITKMASKEDFGTOPOLOGICALNUMBERIMAGEFILTER_H */
