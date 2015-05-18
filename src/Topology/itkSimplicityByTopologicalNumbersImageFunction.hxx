/*=========================================================================
 *
 * This code has been downloaded from 
 * http://www.insight-journal.org/browse/publication/120.
 * 
 * The code was written and contributed by Lamy Julien, 2006.
 * 
 *=========================================================================*/
#ifndef itkSimplicityByTopologicalNumbersImageFunction_hxx
#define itkSimplicityByTopologicalNumbersImageFunction_hxx

#include "itkSimplicityByTopologicalNumbersImageFunction.h"

namespace itk
{

template<typename TImage, typename TForegroundConnectivity, 
         typename TBackgroundConnectivity >
SimplicityByTopologicalNumbersImageFunction<TImage, TForegroundConnectivity, 
                                            TBackgroundConnectivity>
::SimplicityByTopologicalNumbersImageFunction()
  {
  m_TnCounter = TopologicalNumberImageFunction<TImage, 
                  TForegroundConnectivity>::New();
  }

template<typename TImage, typename TForegroundConnectivity, 
        typename TBackgroundConnectivity >
bool
SimplicityByTopologicalNumbersImageFunction<TImage, TForegroundConnectivity, 
                                            TBackgroundConnectivity>
::Evaluate(PointType const & point) const
  {
  typename TImage::IndexType index;
  ConvertPointToNearestIndex(point, index);
  return EvaluateAtIndex(index);
  }


template<typename TImage, typename TForegroundConnectivity, 
         typename TBackgroundConnectivity >
bool
SimplicityByTopologicalNumbersImageFunction<TImage, TForegroundConnectivity, 
                                            TBackgroundConnectivity>
::EvaluateAtIndex(IndexType const & index) const
  {
  std::pair<unsigned char, unsigned char> const result = 
    m_TnCounter->EvaluateAtIndex(index);
  return (result.first==1 && result.second==1);
  }


template<typename TImage, typename TForegroundConnectivity, 
         typename TBackgroundConnectivity >
bool
SimplicityByTopologicalNumbersImageFunction<TImage, TForegroundConnectivity, 
                                            TBackgroundConnectivity>
::EvaluateAtContinuousIndex(ContinuousIndexType const & contIndex) const
  {
  typename TImage::IndexType index;
  ConvertContinuousIndexToNearestIndex(contIndex, index);
  return EvaluateAtIndex(index);
  }

}

#endif // itkSimplicityByTopologicalNumbersImageFunction_hxx
