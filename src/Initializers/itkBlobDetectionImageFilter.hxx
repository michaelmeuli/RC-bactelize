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
#ifndef __ITKBLOBDETECTIONIMAGEFILTER_CX
#define __ITKBLOBDETECTIONIMAGEFILTER_CXX

#include "itkBlobDetectionImageFilter.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
BlobDetectionImageFilter<TInputImage, TOutputImage>
::BlobDetectionImageFilter()
{
    m_MaxDiameter = 20;
    m_MinDiameter = 5;
    m_NumberOfScales = 5;
    
    m_DilateFilter = DilateImageFilterType::New();
}

template <class TInputImage, class TOutputImage>
void
BlobDetectionImageFilter<TInputImage, TOutputImage>::
GenerateData()
{ 

  // Allocate the scale space image (DIM+1)
  typename RealCoDimImageType::Pointer vScaleSpace = RealCoDimImageType::New();
  typename RealCoDimImageType::SizeType vSize;
  typename RealCoDimImageType::RegionType vRegion;
  typename RealImageType::RegionType vSliceRegion = this->GetInput()->GetRequestedRegion();
  for(unsigned int vD =0; vD < TInputImage::ImageDimension; vD++) {
    vSize[vD] = vSliceRegion.GetSize()[vD];
  }
  vSize[TInputImage::ImageDimension] = m_NumberOfScales; // set the last
  vRegion.SetSize(vSize);
  vScaleSpace->SetRegions(vRegion);
  vScaleSpace->Allocate();

  typename LaplacianFilterType::Pointer vLaplacianFilter = LaplacianFilterType::New();
  vLaplacianFilter->SetInput(this->GetInput());
  vLaplacianFilter->SetNormalizeAcrossScale(true);

    itk::ImageRegionIterator<RealCoDimImageType> vCopyIt(vScaleSpace, vScaleSpace->GetLargestPossibleRegion());


    // Generate the scale space image
  float vSigmaMax = 0.5f*m_MaxDiameter/1.414f;
  float vSigmaMin = 0.5f*m_MinDiameter/1.414f;
  float vSigmaInc = (vSigmaMax-vSigmaMin)/(m_NumberOfScales-1);
  float vSigma = vSigmaMin;
  for(unsigned int vI = 0; vI < m_NumberOfScales; vI++) {
    std::cout << "Sigma:\t" <<  vSigma << std::endl;

    vLaplacianFilter->SetSigma(vSigma);
    vLaplacianFilter->Update();

    // copy the slice to the scale space image
    itk::ImageRegionConstIterator<RealImageType> vOrigIt(vLaplacianFilter->GetOutput(), 
            vLaplacianFilter->GetOutput()->GetBufferedRegion());

    for(vOrigIt.GoToBegin(); !vOrigIt.IsAtEnd(); ++vOrigIt) {
      vCopyIt.Set(vOrigIt.Get());
      ++vCopyIt;
    }
    vSigma += vSigmaInc;

  }

  // First do a projection, its faster and the thresholding does not 
  // depend on the number of scale stages.
  typename MinProjIFType::Pointer vMinProjFilter = MinProjIFType::New();
  vMinProjFilter->SetInput(vScaleSpace);
  vMinProjFilter->Update();

  //Find the threshold to select only meaningful/strong local minima
  typename HistogramGeneratorType::Pointer vHistogramGenerator = HistogramGeneratorType::New();
  vHistogramGenerator->SetInput(vMinProjFilter->GetOutput());
  unsigned int vNumberOfBins = 256;
  vHistogramGenerator->SetNumberOfBins(vNumberOfBins);
  vHistogramGenerator->Compute();
    
  typename OtsuThsCalcType::Pointer vOtsuThsCalculator = OtsuThsCalcType::New();
  vOtsuThsCalculator->SetInputHistogram(vHistogramGenerator->GetOutput());
  unsigned int vNumberOfThresholds = 2;
  vOtsuThsCalculator->SetNumberOfThresholds(vNumberOfThresholds);
  vOtsuThsCalculator->Update();
  typename OtsuThsCalcType::OutputType vThs = vOtsuThsCalculator->GetOutput();
    
  typename ThresholdImageFilterType::Pointer vThsIF = ThresholdImageFilterType::New();
  vThsIF->SetInput(vMinProjFilter->GetOutput());
  vThsIF->SetUpper(vThs[vNumberOfThresholds-1]);

  // Within the thresholded image, find local minima
  typename RegMinImageFilterType::Pointer vRegMinFilter = RegMinImageFilterType::New();
  vRegMinFilter->SetInput(vThsIF->GetOutput());
  vRegMinFilter->SetFullyConnected(true);
  vRegMinFilter->SetForegroundValue(1);

  // Create a small blob for each local minimum.
  m_DilateFilter->SetDilateValue(1);
  
  typename SEType::RadiusType vSERadius;
  // TODO: the radius of each blob can be found in the scale space image.
  vSERadius.Fill(m_MinDiameter/3.0f); // such that discs of different blobs do not touch.
  m_DilateFilter->SetKernel(SEType::Ball(vSERadius));
  m_DilateFilter->SetInput(vRegMinFilter->GetOutput());

  // graft the output.
  m_DilateFilter->GraftOutput( this->GetOutput() );
  m_DilateFilter->Update();
  this->GraftOutput(m_DilateFilter->GetOutput() );
}

template <class TInputImage, class TOutputImage>
void
BlobDetectionImageFilter<TInputImage, TOutputImage>::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "Number of scales:" << this->m_NumberOfScales << std::endl;
  os << indent << "Maximal scales:" << this->m_NumberOfScales << std::endl;
  os << indent << "Minimal scale:" << this->m_NumberOfScales << std::endl;
}

} /* end namespace itk */

#endif
