/*=========================================================================
 *
 * This code has been downloaded from 
 * http://www.insight-journal.org/browse/publication/120.
 * 
 * The code was written and contributed by Lamy Julien, 2006.
 * 
 *=========================================================================*/
#ifndef itkTopologicalNumberImageFunction_hxx
#define itkTopologicalNumberImageFunction_hxx

#include "itkTopologicalNumberImageFunction.h"
#include <itkNumericTraits.h>
#include <set>


namespace itk {

    template<typename TImage, typename TFGConnectivity, typename TBGConnectivity >
    TopologicalNumberImageFunction<TImage, TFGConnectivity, TBGConnectivity>
    ::TopologicalNumberImageFunction()
    : m_ComputeForegroundTN(true), m_ComputeBackgroundTN(true) {

        m_IgnoreLabel = itk::NumericTraits<typename TImage::PixelType>::max();
        unsigned int const vImageSize =
                TFGConnectivity::GetInstance().GetNeighborhoodSize();

        m_SubImage = new char[vImageSize];
        m_Offsets = new OffsetType[vImageSize];
        m_DataSubImage = new PixelType[vImageSize];

        // Get the sub-image
        for (unsigned int i = 0; i < vImageSize; ++i) {
            int remainder = i;

            // Get current offset
            for (unsigned int j = 0; j < TFGConnectivity::Dimension; ++j) {
                m_Offsets[i][j] = remainder % 3;
                remainder -= m_Offsets[i][j];
                remainder /= 3;
                --m_Offsets[i][j];
            }
        }
    }

    template<typename TImage, typename TFGConnectivity, typename TBGConnectivity >
    TopologicalNumberImageFunction<TImage, TFGConnectivity, TBGConnectivity>
    ::~TopologicalNumberImageFunction() {
       delete m_SubImage;
       delete m_Offsets;
       delete m_DataSubImage;
    }

    template<typename TImage, typename TFGConnectivity, typename TBGConnectivity >
    std::pair<unsigned int, unsigned int>
    TopologicalNumberImageFunction<TImage, TFGConnectivity, TBGConnectivity>
    ::Evaluate(PointType const & point) const {
        typename TImage::IndexType index;
        this->ConvertPointToNearestIndex(point, index);
        return EvaluateAtIndex(index);
    }

    
    /// Evaluates the topological numbers for the current label at index. 
    template<typename TImage, typename TFGConnectivity, typename TBGConnectivity >
    std::pair<unsigned int, unsigned int>  
    TopologicalNumberImageFunction<TImage, TFGConnectivity, TBGConnectivity>
    ::EvaluateFGTNAtIndex(IndexType const & index) const {
        return EvaluateFGTNOfLabelAtIndex(index, this->GetInputImage()->GetPixel(index));
    }
    
    /// Evaluates the topological number for a particular label at the index 
    /// provided.
    template<typename TImage, typename TFGConnectivity, typename TBGConnectivity >
    std::pair<unsigned int, unsigned int>  
    TopologicalNumberImageFunction<TImage, TFGConnectivity, TBGConnectivity>
    ::EvaluateFGTNOfLabelAtIndex(IndexType const & index, PixelType aLabel) const {
        unsigned int const imageSize =
        TFGConnectivity::GetInstance().GetNeighborhoodSize();
        
        readDataSubImage(index);
        
        for (unsigned int i = 0; i < imageSize; ++i) {
            m_SubImage[i] = (m_DataSubImage[i] == aLabel) ? 255 : 0;
        }
        
        unsigned int const middle = imageSize / 2;
        m_SubImage[middle] = 0;
        
        // Topological number in the foreground
        FGandBGTopoNbPairType vFGBGTopoPair;
        m_ForegroundUnitCubeCCCounter.SetImage(m_SubImage, m_SubImage + imageSize);
        
        vFGBGTopoPair.first = m_ForegroundUnitCubeCCCounter();
        
        // Invert the sub-image
        for (unsigned int bit = 0; bit < middle; ++bit) {
            m_SubImage[bit] = 255 - m_SubImage[bit];
        }
        for (unsigned int bit = middle; bit < imageSize - 1; ++bit) {
            m_SubImage[bit + 1] = 255 - m_SubImage[bit + 1];
        }
        
        // Topological number in the background
        m_BackgroundUnitCubeCCCounter.SetImage(m_SubImage,
                m_SubImage + TBGConnectivity::GetInstance().GetNeighborhoodSize());
        
        vFGBGTopoPair.second = m_BackgroundUnitCubeCCCounter();
        
        return vFGBGTopoPair;
    }
    
    /// Evaluates the topological number for all adjacent (to index) labels when
    /// considering all other labels as background. 
    template<typename TImage, typename TFGConnectivity, typename TBGConnectivity >
    std::vector<std::pair<unsigned int, std::pair<unsigned int, unsigned int> > >
    TopologicalNumberImageFunction<TImage, TFGConnectivity, TBGConnectivity>
    ::EvaluateAdjacentRegionsFGTNAtIndex(IndexType const & index) const {

        ForegroundTopologicalNumbersType vTNvector;

        unsigned int const imageSize =
                TFGConnectivity::GetInstance().GetNeighborhoodSize();

        typedef std::set<PixelType> AdjacentLabelsSetType;
        AdjacentLabelsSetType vAdjacentLabels;

        readDataSubImage(index);
        
        static const OffsetType * const vFGNeighbors =
                TFGConnectivity::GetInstance().GetNeighborsITKOffsets();
        static const unsigned int vFGNeighborhoodsize =
                TFGConnectivity::GetInstance().GetNumberOfNeighbors();

        for (unsigned int vN = 0; vN < vFGNeighborhoodsize; vN++) {
            int vLinearIndex = 0;
            int vBase3 = 1;
            for (unsigned int vDim = 0; vDim < TFGConnectivity::Dimension; ++vDim) {
                vLinearIndex += (vFGNeighbors[vN][vDim]+1) * vBase3;
                vBase3 *= 3;
            }

            if (m_DataSubImage[vLinearIndex] != itk::NumericTraits<typename TImage::PixelType>::Zero) {
                vAdjacentLabels.insert(m_DataSubImage[vLinearIndex]);
            }
        }

        typedef typename AdjacentLabelsSetType::iterator AdjacentLabelsSetIteratorType;
        AdjacentLabelsSetIteratorType vLabelsIt = vAdjacentLabels.begin();
        AdjacentLabelsSetIteratorType vLabelsEnd = vAdjacentLabels.end();

        for (; vLabelsIt != vLabelsEnd; ++vLabelsIt) {
            FGandBGTopoNbPairType vPair = EvaluateFGTNOfLabelAtIndex(index, *vLabelsIt);
            vTNvector.push_back(ForegroundTopologicalNumberType(*vLabelsIt, vPair));
        }
        
        return vTNvector;
    }

    /// Classical DT topological number
    template<typename TImage, typename TFGConnectivity, typename TBGConnectivity >
    std::pair<unsigned int, unsigned int>
    TopologicalNumberImageFunction<TImage, TFGConnectivity, TBGConnectivity>
    ::EvaluateAtIndex(IndexType const & index) const {
        unsigned int const imageSize =
                TFGConnectivity::GetInstance().GetNeighborhoodSize();

        
        readDataSubImage(index);

        // Get the binary sub-image
        for (unsigned int i = 0; i < imageSize; ++i) {

            typename TImage::PixelType vImgVal = m_DataSubImage[i];
            m_SubImage[i] =
                    (vImgVal != itk::NumericTraits<typename TImage::PixelType>::Zero) ?
                    255 : 0;

        }

        unsigned int const middle = imageSize / 2;

        m_SubImage[middle] = 0;

        // Topological number in the foreground
        m_ForegroundUnitCubeCCCounter.SetImage(m_SubImage, m_SubImage + imageSize);
        unsigned int const ccNumber =
                m_ComputeForegroundTN ? m_ForegroundUnitCubeCCCounter() : 0;

        // Invert the sub-image
        for (unsigned int bit = 0; bit < middle; ++bit) {
            m_SubImage[bit] = 255 - m_SubImage[bit];
        }
        for (unsigned int bit = middle; bit < imageSize - 1; ++bit) {
            m_SubImage[bit + 1] = 255 - m_SubImage[bit + 1];
        }

        // Topological number in the background
        m_BackgroundUnitCubeCCCounter.SetImage(m_SubImage,
                m_SubImage + TBGConnectivity::GetInstance().GetNeighborhoodSize());
        assert(TFGConnectivity::GetInstance().GetNeighborsPoints());

        unsigned int const backgroundCcNumber =
                m_ComputeBackgroundTN ? m_BackgroundUnitCubeCCCounter() : 0;        

        return std::pair<unsigned int, unsigned int>(ccNumber, backgroundCcNumber);
    }

    template<typename TImage, typename TFGConnectivity, typename TBGConnectivity >
    std::pair<unsigned int, unsigned int>
    TopologicalNumberImageFunction<TImage, TFGConnectivity, TBGConnectivity>
    ::EvaluateAtContinuousIndex(ContinuousIndexType const & contIndex) const {
        typename TImage::IndexType index;
        this->ConvertContinuousIndexToNearestIndex(contIndex, index);
        return EvaluateAtIndex(index);
    }

    template<typename TImage, typename TFGConnectivity, typename TBGConnectivity >
    void
    TopologicalNumberImageFunction<TImage, TFGConnectivity, TBGConnectivity>
    ::readDataSubImage(IndexType const & index) const{
        unsigned int const imageSize =
                TFGConnectivity::GetInstance().GetNeighborhoodSize();

        for (unsigned int i = 0; i < imageSize; ++i) {
            m_DataSubImage[i] = abs(this->GetInputImage()->GetPixel(index + m_Offsets[i]));
            if (m_DataSubImage[i] == m_IgnoreLabel) m_DataSubImage[i] =
                    itk::NumericTraits<typename TImage::PixelType>::Zero;
        }
        
    }
    
    template<typename TImage, typename TFGConnectivity, typename TBGConnectivity >
            UnitCubeCCCounter< TFGConnectivity >
            TopologicalNumberImageFunction<TImage, TFGConnectivity, TBGConnectivity>
            ::m_ForegroundUnitCubeCCCounter = UnitCubeCCCounter< TFGConnectivity >();


    template<typename TImage, typename TFGConnectivity, typename TBGConnectivity >
            UnitCubeCCCounter< TBGConnectivity >
            TopologicalNumberImageFunction<TImage, TFGConnectivity, TBGConnectivity>
            ::m_BackgroundUnitCubeCCCounter = UnitCubeCCCounter< TBGConnectivity >();

}

#endif // itkTopologicalNumberImageFunction_hxx
