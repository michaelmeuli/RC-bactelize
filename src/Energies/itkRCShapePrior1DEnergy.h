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
#ifndef ITKRCSHAPEPRIOR1DENERGY_H
#define	ITKRCSHAPEPRIOR1DENERGY_H

#include "itkRCInternalEnergyBaseClass.h"
#include <itkNumericTraits.h>

#include "itkOrderedPoint.h"
#include "itkMoments1D.h"
#include "itkMoments1DEnergy.h"
#include "itkFixedArray.h"
#include "itkImageRegionConstIteratorWithIndex.h"

namespace itk
{

/**
 * @brief Returns a regularizing energy difference based L1 distances of moment
 * vectors of an object and a reference object of a certain order.
 * 
 */
template<typename TLabelImage, typename TTemplateImage, typename TEnergyDifference = float>
class ITK_EXPORT RCShapePrior1DEnergy :
  public RCInternalEnergyBaseClass<TLabelImage, TEnergyDifference>
  {
  public :
    /*
     * Standard ITK declarations
     */
    typedef RCShapePrior1DEnergy<TLabelImage, TTemplateImage, TEnergyDifference> Self;
    typedef RCInternalEnergyBaseClass<TLabelImage, TEnergyDifference> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<Self const> ConstPointer;

    itkNewMacro(Self);
    itkTypeMacro(RCShapePrior1DEnergy, RCInternalEnergyBaseClass);
    
    typedef TTemplateImage TemplateImageType;
    typedef typename TemplateImageType::PixelType TemplatePixelType;
    typedef typename TemplateImageType::Pointer TemplateImagePointerType;
    typedef typename TemplateImageType::ConstPointer TemplateImageConstPointerType;
    
    typedef typename Superclass::IndexType IndexType;
    typedef typename Superclass::LabelPixelType LabelPixelType;
    typedef typename Superclass::LabelAbsPixelType LabelAbsPixelType;
    typedef typename Superclass::EnergyDifferenceType EnergyDifferenceType;
    typedef typename Superclass::InternalEnergyReturnType InternalEnergyReturnType;
    typedef typename Superclass::LabelImageType LabelImageType;
    typedef typename Superclass::OffsetType OffsetType;
    typedef typename Superclass::SizeType SizeType;
    typedef typename Superclass::RegionType RegionType;
    
    typedef FixedArray<float, LabelImageType::ImageDimension+1> Centroid_Coeff_Type;
    
    inline InternalEnergyReturnType EvaluateEnergyDifference(IndexType aInd, 
            LabelPixelType aLabelBefore, 
            LabelPixelType aLabelAfter);
   
    /** Add a point to the statistics of region. */
    inline void AddPoint(IndexType aInd, LabelAbsPixelType aRegion);
    
    /** Remove a point from the statistics of a region. */
    inline void RemovePoint(IndexType aInd /* = 0 */, LabelAbsPixelType aRegion);
    
    /** Implement the virtual method to clean-up the region statistics. */
    void KillRegion(LabelAbsPixelType aL);

    /** 
     * Sets the value which is ignored in the template image.
     * All the values except the background value are considered foreground and
     * hence belong to the template object. All foreground pixel contribute to
     * the moment statistics. Defaults to 0.
     * See also SetBackgroundLabel
     */
    itkSetMacro(BackgroundValue, TemplatePixelType);
    /** Get the pixel value that is considered background in the template image. */
    itkGetConstMacro(BackgroundValue, TemplatePixelType);
    
    /** Set the label value which is considered to be background. For this label
     no shape prior will be imposed. Usually the background label is zero.*/
    itkSetMacro(BackgroundLabel, LabelAbsPixelType);
    
    /** Get the background label for which no shape prior will be imposed. */
    itkGetMacro(BackgroundLabel, LabelAbsPixelType);

    /** Set the order of the moments considered to measure moment distances. 
     20 to 40 are common values. Defaults to 20. */
    itkSetMacro(OrderOfMoments, unsigned int);

    /** Get the order of the moments considered to measure moment distances. */
    itkGetConstMacro(OrderOfMoments, unsigned int);

    /** 
     * Sets how many moment distributions should be considered. Each of them
     * has a different reference point. The first reference point is positioned
     * at the center of mass. The following are positioned at in direction of
     * the principal components. 
     */
    itkSetMacro(ShapePriorNbCentroids, unsigned int);
    
    /**
     * Sets how many moment distributions should be considered. Each of them
     * having a difference reference point.
     */
    itkGetMacro(ShapePriorNbCentroids, unsigned int);
    
    /** Set the coefficients for the different moment vectors. Each moment 
     * vector is generated w.r.t. different reference positions. With 
     * SetShapePriorCentroidCoeff the influence of each might be changed. 
     */
    itkSetMacro(ShapePriorCentroidCoeff,Centroid_Coeff_Type);
    
    /** Get the coefficients for the different moment vectors. Each moment 
     * vector is generated w.r.t. different reference positions.
     */
    itkGetConstMacro(ShapePriorCentroidCoeff,Centroid_Coeff_Type);
    
    /** Set the minimum size of the region for which the shape prior kicks in. 
     * For too small regions the discretization noise is too large and the prior
     * might work badly.
     */
    itkSetMacro(ShapePriorMinSize, unsigned int);
    
    /** Get the minimum area for which the prior kicks in. */
    itkGetMacro(ShapePriorMinSize, unsigned int);
    
    void PrepareEnergyCalculation();
    void PrepareEnergyCalculationForIteration();
    
    /** Get the shape template image. */
    const TemplateImageType * GetTemplateImage() const
    { return m_TemplateImage.GetPointer(); }
    
    /** Set the shape template image. This needs to be called to define a valid 
     *  energy.
     */
    void SetTemplateImage(const TemplateImageType *aPtr)
    { m_TemplateImage = aPtr; }
  
    
    protected:
        RCShapePrior1DEnergy();
        ~RCShapePrior1DEnergy() {};
        void PrintSelf(std::ostream & os, Indent indent) const;
        
        
    private:
        RCShapePrior1DEnergy(const Self &);  //purposely not implemented
        void operator=(const Self &); //purposely not implemented
                    
        TemplatePixelType m_BackgroundValue;
        LabelAbsPixelType m_BackgroundLabel;
        TemplateImageConstPointerType m_TemplateImage;        
        unsigned int m_OrderOfMoments;
        
        typedef OrderedPoint<double, TemplateImageType::ImageDimension> MomentsPointType;
        typedef Moments1DEnergy<MomentsPointType> Moments1DType;
        typedef typename Moments1DType::Pointer Moments1DPointerType;
        Moments1DPointerType m_Moments1DReference;

        typedef std::map<LabelAbsPixelType, Moments1DPointerType> Moments1DStatisticsType;
        typedef typename Moments1DStatisticsType::iterator Moments1DStatisticsIteratorType;
        Moments1DStatisticsType m_Moments1DStatistics;
        
        unsigned int m_ShapePriorNbCentroids;
        Centroid_Coeff_Type m_ShapePriorCentroidCoeff;
        unsigned int m_ShapePriorMinSize;
        
}; // end class
} // end namespace


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRCShapePrior1DEnergy.hxx"
#endif
 
#endif	/* ITKRCSHAPEPRIOR1DENERGY_H */

