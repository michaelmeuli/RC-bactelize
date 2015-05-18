/*=========================================================================
 *
 * This code has been downloaded from 
 * http://www.insight-journal.org/browse/publication/120.
 * 
 * The code was written and contributed by Lamy Julien, 2006.
 * 
 *=========================================================================*/
#ifndef itkTopologicalNumberImageFunction_h
#define itkTopologicalNumberImageFunction_h


#include <utility>
#include <vector>

#include <itkImageFunction.h>

#include "itkBackgroundConnectivity.h"
#include "itkUnitCubeCCCounter.h"


namespace itk {

    /**
     * @brief Compute the topological numbers of an image at given index.
     *
     * Topological numbers characterize the topological properties of a point. They
     * are defined in an article by G. Bertrand and G. Malandain : "A new
     * characterization of three-dimensional simple points"; Pattern Recognition
     * Letters; 15:169--175; 1994.
     */
    template<typename TImage,
            typename TFGConnectivity,
            typename TBGConnectivity =
            typename BackgroundConnectivity<TFGConnectivity>::Type >
            class ITK_EXPORT TopologicalNumberImageFunction :
            public itk::ImageFunction<TImage, std::pair<unsigned int, unsigned int> > {
    public:
        /**
         * @name Standard ITK declarations
         */
        //@{
        typedef TopologicalNumberImageFunction<TImage, TFGConnectivity> Self;
        typedef itk::ImageFunction<TImage, std::pair<unsigned int, unsigned int> >
        Superclass;
        typedef itk::SmartPointer<Self> Pointer;
        typedef itk::SmartPointer<Self const> ConstPointer;

        typedef itk::Offset<TFGConnectivity::Dimension> OffsetType;
        typedef typename TImage::PixelType PixelType;

        itkNewMacro(Self);
        itkTypeMacro(TopologicalNumberImageFunction, ImageFunction);

        typedef typename Superclass::PointType PointType;
        typedef typename Superclass::ContinuousIndexType ContinuousIndexType;
        typedef typename Superclass::IndexType IndexType;
        //@}
        
        typedef std::pair<unsigned int, unsigned int> FGandBGTopoNbPairType;
        typedef std::pair<unsigned int, FGandBGTopoNbPairType> ForegroundTopologicalNumberType;        
        typedef std::vector<ForegroundTopologicalNumberType> ForegroundTopologicalNumbersType;


        /**
         * @brief Initialize the functor so that the topological numbers are
         * computed for both the foreground and the background.
         */
        TopologicalNumberImageFunction();
        ~TopologicalNumberImageFunction();

        /**
         * @name Evaluation functions
         *
         * These functions evaluate the topological number at the index.
         */
        //@{
        FGandBGTopoNbPairType
        Evaluate(PointType const & point) const;

        FGandBGTopoNbPairType
        EvaluateAtIndex(IndexType const & index) const;

        FGandBGTopoNbPairType
        EvaluateAtContinuousIndex(ContinuousIndexType const & contIndex) const;
     
        FGandBGTopoNbPairType 
        EvaluateFGTNAtIndex(IndexType const & index) const;
        
        FGandBGTopoNbPairType 
        EvaluateFGTNOfLabelAtIndex(IndexType const & index, PixelType aLabel) const;

        ForegroundTopologicalNumbersType
        EvaluateAdjacentRegionsFGTNAtIndex(IndexType const & index) const;

    private:
        void readDataSubImage(IndexType const & index) const;
        
    public:
        //@}

        /**
         * @name Selectors for the computation of fore- and background topological
         * numbers.
         *
         * These two members allow to selectively compute the topological
         * numbers for the background and the foreground. They are both set to true
         * during the construction of the object.
         */
        //@{
        itkGetConstMacro(ComputeForegroundTN, bool);
        itkSetMacro(ComputeForegroundTN, bool);
        itkGetConstMacro(ComputeBackgroundTN, bool);
        itkSetMacro(ComputeBackgroundTN, bool);
        //@}

    private:
        static UnitCubeCCCounter< TFGConnectivity > m_ForegroundUnitCubeCCCounter;
        static UnitCubeCCCounter< TBGConnectivity > m_BackgroundUnitCubeCCCounter;

        TopologicalNumberImageFunction(Self const &); // not implemented
        Self & operator=(Self const &); // not implemented

        bool m_ComputeForegroundTN;
        bool m_ComputeBackgroundTN;

        char* m_SubImage;
        OffsetType* m_Offsets;
        PixelType* m_DataSubImage;

        typename TImage::PixelType m_IgnoreLabel;
    };

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTopologicalNumberImageFunction.hxx"
#endif

#endif // itkTopologicalNumberImageFunction_h
