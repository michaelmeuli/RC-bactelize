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

#ifndef ITKRCUTILS_H
#define	ITKRCUTILS_H

#include "itkConnectivity.h"
#include "itkBackgroundConnectivity.h"

namespace itk {
    

/** \class RCUtils
 * \brief A convenience class for the region competition framework.
 * 
 * The first template argument specifies the label image type used, which must
 * be of integral type. The second argument specifies the foreground connectivity
 * type used. The background connectivity type is selected using the digital
 * topology framework by J.Lamy.
 * 
 * \ingroup Segmentation
 */
template< typename TLabelImage, unsigned int VCellDim >
class RCUtils
{
public:
  typedef RCUtils<TLabelImage, VCellDim> Self;
  typedef TLabelImage LabelImageType;
  typedef typename LabelImageType::IndexType IndexType;
  typedef typename LabelImageType::SizeType SizeType;
  typedef typename LabelImageType::RegionType RegionType;
  typedef typename LabelImageType::Pointer LabelImagePointerType;
  typedef typename LabelImageType::IndexType LabelImageIndexType;
  typedef typename LabelImageType::PixelType LabelPixelType;
  typedef typename LabelImageType::OffsetType LabelImageOffsetType;
  typedef typename LabelImageType::OffsetType OffsetType;
  
  
  // TODO: get the integral type of LabelPixelType:
  typedef unsigned int LabelAbsPixelType;
  
  
  static inline bool IsEnclosedByLabel_BGConnectivity(LabelImagePointerType aLabelImage,
          const LabelImageIndexType& aIndex, 
          LabelAbsPixelType aLabel);
  
  static inline bool IsSingleFGPoint(LabelImagePointerType aLabelImage,
          const LabelImageIndexType& aIndex, LabelAbsPixelType aAbsLabel);
  
  static inline bool IsBoundaryPoint(LabelImagePointerType aLabelImage, 
          const LabelImageIndexType& aIndex);
          
  /** Connectivity, Topology types and their members */
  typedef Connectivity<LabelImageType::ImageDimension, LabelImageType::ImageDimension - 1 > ForegroundConnectivityType;
  typedef typename BackgroundConnectivity<ForegroundConnectivityType>::Type BackgroundConnectivityType;
  

  // TODO: the GetNeighborsITKOffsets() return a pointer to an array that we 
  //       have to delete, which we don't -> memory leaks.
  
  static unsigned int m_NeighborhoodSize_FG_Connectivity;
  
  static const OffsetType * m_NeighborsOffsets_FG_Connectivity;
    
  static unsigned int m_NeighborhoodSize_BG_Connectivity;

  static const OffsetType * m_NeighborsOffsets_BG_Connectivity;
  

//  template<ImageType>
//  static void WriteImageToDisk(ImageType::Pointer aImgPtr, const char* aFilename , float aScale);
        
}; // end class

} // end namespace itk

#include "itkRCUtils.hxx"

#endif	/* ITKRCUTILS_H */

