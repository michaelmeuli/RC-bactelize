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
#ifndef ITKRCSHAPEPRIOR23DENERGY_H
#define	ITKRCSHAPEPRIOR23DENERGY_H

#include "itkRCInternalEnergyBaseClass.h"
#include <itkNumericTraits.h>
#include <map>
#include "Moments2D.h"
//#include "Moments3D.h" // ?

#include "itkFixedArray.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "vnl/vnl_matrix.h"

namespace itk
{

/**
 * @brief Returns a regularizing energy difference based L1 distances of moment
 * vectors of an object and a reference object of a certain order.
 * 
 */
template<typename TLabelImage, typename TTemplateImage, typename TEnergyDifference = float>
class ITK_EXPORT RCShapePrior23DEnergy :
  public RCInternalEnergyBaseClass<TLabelImage, TEnergyDifference>
  {
  public :
    /*
     * Standard ITK declarations
     */
    typedef RCShapePrior23DEnergy<TLabelImage, TTemplateImage, TEnergyDifference> Self;
    typedef RCInternalEnergyBaseClass<TLabelImage, TEnergyDifference> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<Self const> ConstPointer;

    itkNewMacro(Self);
    itkTypeMacro(RCShapePrior23DEnergy, RCInternalEnergyBaseClass);
    
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
        RCShapePrior23DEnergy();
        ~RCShapePrior23DEnergy() {};
        void PrintSelf(std::ostream & os, Indent indent) const;
        
        
    private:
        RCShapePrior23DEnergy(const Self &);  //purposely not implemented
        void operator=(const Self &); //purposely not implemented
                    
        void ShapePriorReferenceInit();
        void InitializeShapePrior();
        InternalEnergyReturnType CalcFullEnergy();
        
        TemplatePixelType m_BackgroundValue;
        LabelAbsPixelType m_BackgroundLabel;
        TemplateImageConstPointerType m_TemplateImage;        
        unsigned int m_OrderOfMoments;
        
        typedef Moments2D Moments2DType;
        typedef Moments2DType::Pointer Moments2DPointerType;
        typedef std::map<LabelAbsPixelType, Moments2DPointerType> LabelToMoments2DMapType;
        std::vector<IndexType> m_PointsToBeAdded;
        std::vector<IndexType> m_PointsToBeDeleted;
        LabelToMoments2DMapType m_Moments2Ds;
        
        typedef std::map<LabelAbsPixelType, unsigned int>  CountStatisticsMapType;
        typedef typename CountStatisticsMapType::iterator CountStatisticsIteratorType;
        
        /** The number of pixel within a regin. */
        CountStatisticsMapType m_Count;
        

        
}; // end class
} // end namespace


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRCShapePrior23DEnergy.hxx"
#endif
 
#endif	/* ITKRCSHAPEPRIOR23DENERGY_H */

