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
/* 
 * File:   itkEnergyFunction.h
 * Author: janickc
 *
 * Created on April 9, 2010, 4:33 PM
 */

#ifndef _ITKENERGYFUNCTION_H
#define	_ITKENERGYFUNCTION_H

#include "itkFunctionBase.h"

namespace itk
{


/** \class EnergyFunction
 * \brief Evaluates a function of an image at specified position.
 *
 * ImageFunction is a baseclass for all objects that evaluates
 * a function of an image at an image iterator.
 *
 * The input image is set via method SetInputImage().
 * Evaluate (IteratorType) returns a number of ValueType type.
 *
 */
template <
class TInputImage,
class TOutput = float,
>
class ITK_EXPORT EnergyFunction : public FunctionBase< ImageIterator<TInputImage>,
                       TOutput >
{
public:
  /** Dimension underlying input image. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Standard class typedefs. */
  typedef ImageFunction                                         Self;
  typedef FunctionBase<ImageIterator<TInputImage>,TOutput       Superclass;
  typedef SmartPointer<Self>                                    Pointer;
  typedef SmartPointer<const Self>                              ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(EnergyFunction, FunctionBase);

  /** InputImageType typedef support. */
  typedef TInputImage InputImageType;

  /** InputPixel typedef support */
  typedef typename InputImageType::PixelType InputPixelType;

  /** InputImagePointer typedef support */
  typedef typename InputImageType::ConstPointer InputImageConstPointer;

  /** OutputType typedef support. */
  typedef TOutput OutputType;

  /** CoordRepType typedef support. */
  typedef TCoordRep CoordRepType;

  /** Index Type. */
  typedef typename InputImageType::IndexType IndexType;


  /** Set the input image.
   * \warning this method caches BufferedRegion information.
   * If the BufferedRegion has changed, user must call
   * SetInputImage again to update cached values. */
  virtual void SetInputImage( const InputImageType * ptr );

  /** Get the input image. */
  const InputImageType * GetInputImage() const
    { return m_Image.GetPointer(); }

  /** Evaluate the function at the current iterator position */
  virtual TOutput Evaluate( const ImageIterator& aIterator ) const = 0;

  /** Evaluate the function at specified Index position.
   * Subclasses must provide this method. */
  virtual TOutput EvaluateAtIndex( const IndexType & index ) const = 0;

  /** Evaluate the function at specified ContinuousIndex position.
   * Subclasses must provide this method. */
  virtual TOutput EvaluateAtContinuousIndex(
    const ContinuousIndexType & index ) const = 0;

  /** Check if an index is inside the image buffer.
   * If ITK_USE_CENTERED_PIXEL_COORDINATES_CONSISTENTLY is on,
   * we take into account the fact that each voxel has its
   * center at the integer coordinate and extends half way
   * to the next integer coordinate.
   * \warning For efficiency, no validity checking of
   * the input image is done. */
  virtual bool IsInsideBuffer( const IndexType & index ) const
    {
    for( unsigned int j = 0; j < ImageDimension; j++ )
      {
      if( index[j] < m_StartIndex[j] )
        {
        return false;
        }
      if( index[j] > m_EndIndex[j] )
        {
        return false;
        }
      }
    return true;
    }

  /** Check if a continuous index is inside the image buffer.
   * \warning For efficiency, no validity checking of
   * the input image is done. */
  virtual bool IsInsideBuffer( const ContinuousIndexType & index ) const
    {
    for( unsigned int j = 0; j < ImageDimension; j++ )
      {
      if( index[j] < m_StartContinuousIndex[j] )
        {
        return false;
        }
#ifdef ITK_USE_CENTERED_PIXEL_COORDINATES_CONSISTENTLY
      if( index[j] >= m_EndContinuousIndex[j] )
#else
      if( index[j] > m_EndContinuousIndex[j] )
#endif
      {
        return false;
        }
      }
    return true;
    }

  /** Check if a point is inside the image buffer.
   * \warning For efficiency, no validity checking of
   * the input image pointer is done. */
  virtual bool IsInsideBuffer( const PointType & point ) const
    {
    ContinuousIndexType index;
    m_Image->TransformPhysicalPointToContinuousIndex( point, index );
    return this->IsInsideBuffer( index );
    }

  /** Convert point to nearest index. */
  void ConvertPointToNearestIndex( const PointType & point,
    IndexType & index ) const
    {
    ContinuousIndexType cindex;
    m_Image->TransformPhysicalPointToContinuousIndex( point, cindex );
    this->ConvertContinuousIndexToNearestIndex( cindex, index );
    }

  /** Convert point to continuous index */
  void ConvertPointToContinousIndex( const PointType & point,
    ContinuousIndexType & cindex ) const
    {
    itkWarningMacro("Please change your code to use ConvertPointToContinuousIndex "
      << "rather than ConvertPointToContinousIndex. The latter method name was "
      << "mispelled and the ITK developers failed to correct it before it was released."
      << "The mispelled method name is retained in order to maintain backward compatibility.");
    this->ConvertPointToContinuousIndex( point, cindex );
    }

  /** Convert point to continuous index */
  void ConvertPointToContinuousIndex( const PointType & point,
    ContinuousIndexType & cindex ) const
    {
    m_Image->TransformPhysicalPointToContinuousIndex( point, cindex );
    }

  /** Convert continuous index to nearest index. */
  inline void ConvertContinuousIndexToNearestIndex( const ContinuousIndexType & cindex,
    IndexType & index ) const
    {
    index.CopyWithRound( cindex );
    }

  itkGetConstReferenceMacro(StartIndex, IndexType);
  itkGetConstReferenceMacro(EndIndex, IndexType);

  itkGetConstReferenceMacro(StartContinuousIndex, ContinuousIndexType);
  itkGetConstReferenceMacro(EndContinuousIndex, ContinuousIndexType);

protected:
  ImageFunction();
  ~ImageFunction() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Const pointer to the input image. */
  InputImageConstPointer  m_Image;

  /** Cache some values for testing if indices are inside buffered region. */
  IndexType               m_StartIndex;
  IndexType               m_EndIndex;
  ContinuousIndexType     m_StartContinuousIndex;
  ContinuousIndexType     m_EndContinuousIndex;

private:
  ImageFunction(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

}// end namespace itk


// Define instantiation macro for this template.
#define ITK_TEMPLATE_ImageFunction(_, EXPORT, x, y) namespace itk { \
  _(3(class EXPORT ImageFunction< ITK_TEMPLATE_3 x >)) \
  namespace Templates { typedef ImageFunction< ITK_TEMPLATE_3 x > ImageFunction##y; } \
  }

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkImageFunction+-.h"
#endif

#if ITK_TEMPLATE_HXX
//# include "itkEnergyFunction.hxx"
#endif

#endif	/* _ITKENERGYFUNCTION_H */

