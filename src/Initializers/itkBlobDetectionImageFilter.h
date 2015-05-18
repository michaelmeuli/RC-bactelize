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
#ifndef __ITKBLOBDETECTIONIMAGEFILTER_H
#define __ITKBLOBDETECTIONIMAGEFILTER_H

#include "itkAddImageFilter.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkRegionalMinimaImageFilter.h" 
#include "itkAndImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkMinimumProjectionImageFilter.h"
#include "itkHistogram.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkOtsuMultipleThresholdsCalculator.h"

namespace itk {

template <class TInputImage, class TOutputImage>
class ITK_EXPORT BlobDetectionImageFilter :
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public:

  typedef BlobDetectionImageFilter               Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage> Superclass;
  typedef SmartPointer<Self>                        Pointer;
  typedef SmartPointer<const Self>                  ConstPointer;


  /** Method for creation through object factory */
  itkNewMacro(Self);

  /** Run-time type information */
  itkTypeMacro(BlobDetectionImageFilter, ImageToImageFilter);

  /** Display */
  void PrintSelf( std::ostream& os, Indent indent ) const;

  typedef typename TInputImage::PixelType InputPixelType;
  typedef typename TOutputImage::PixelType OutputPixelType;

  typedef float RealPixelType;
  typedef itk::Image<RealPixelType, TInputImage::ImageDimension> RealImageType;
  typedef itk::Image<RealPixelType, TInputImage::ImageDimension+1> RealCoDimImageType;  
  typedef unsigned int UIntPixelType;
  typedef itk::Image<UIntPixelType, TInputImage::ImageDimension> UIntImageType;
  

  itkGetMacro(NumberOfScales, unsigned int);
  itkSetMacro(NumberOfScales, unsigned int);
  itkGetMacro(MinDiameter, float);
  itkSetMacro(MinDiameter, float);
  itkGetMacro(MaxDiameter, float);
  itkSetMacro(MaxDiameter, float);

  
protected:

  BlobDetectionImageFilter();

protected:

    typedef itk::LaplacianRecursiveGaussianImageFilter<TInputImage, RealImageType> 
    LaplacianFilterType;
    
    typedef itk::Statistics::ScalarImageToHistogramGenerator<RealImageType>
    HistogramGeneratorType;
    typedef typename HistogramGeneratorType::HistogramType HistogramType;
    typedef  itk::OtsuMultipleThresholdsCalculator<HistogramType> OtsuThsCalcType;
    
    typedef itk::MinimumProjectionImageFilter<RealCoDimImageType, RealImageType>
    MinProjIFType;

    typedef itk::ThresholdImageFilter<RealImageType> ThresholdImageFilterType;
    
    typedef itk::FlatStructuringElement<TInputImage::ImageDimension> SEType;
    
    typedef itk::RegionalMinimaImageFilter<RealImageType, UIntImageType> RegMinImageFilterType;
    
    typedef itk::BinaryDilateImageFilter<UIntImageType, UIntImageType, SEType>
    DilateImageFilterType;
    
    typedef itk::AndImageFilter<UIntImageType, UIntImageType> AndFilterType;
    
    typedef itk::AddImageFilter<UIntImageType, UIntImageType> 
    AddFilterType;

    typedef itk::BinaryThresholdImageFilter<UIntImageType, TOutputImage> BinaryThsFilterType;
    
    void GenerateData();

private:

  BlobDetectionImageFilter(Self&);   // intentionally not implemented
  void operator=(const Self&);          // intentionally not implemented

  unsigned int m_NumberOfScales;
  float m_MaxDiameter;
  float m_MinDiameter;
  
//  typename UIntImageType::Pointer m_SumOfConsecDetections;
  
  typename DilateImageFilterType::Pointer m_DilateFilter;
  

};

} /* namespace itk */

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlobDetectionImageFilter.hxx"
#endif

#endif
