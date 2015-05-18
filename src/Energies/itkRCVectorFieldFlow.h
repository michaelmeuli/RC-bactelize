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
#ifndef ITKRCVECTORFIELDFLOW_H
#define	ITKRCVECTORFIELDFLOW_H

#include "itkRCExternalEnergyBaseClass.h"
#include "itkConstNeighborhoodIterator.h"

#include <map>
namespace itk
{

/**
 * @brief Computes energy ...
 */
template<typename TLabelImage, typename TDataImage, typename TVectorImage, typename TEnergyDifference = float>
class ITK_EXPORT RCVectorFieldFlow :
  public RCExternalEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
  {
  public :
    /**
     * @name Standard ITK declarations
     */
    //@{
    typedef RCVectorFieldFlow<TLabelImage, TDataImage, TVectorImage, TEnergyDifference> Self;
    typedef RCExternalEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<Self const> ConstPointer;

    itkNewMacro(Self);
    itkTypeMacro(RCVectorFieldFlow, RCExternalEnergyBaseClass);

    typedef TVectorImage VectorImageType;
    typedef typename VectorImageType::ConstPointer VectorImageConstPointerType;
    
    typedef typename Superclass::IndexType IndexType;
    typedef typename Superclass::LabelPixelType LabelPixelType;
    typedef typename Superclass::LabelAbsPixelType LabelAbsPixelType;
    typedef typename Superclass::DataPixelType DataPixelType;
    typedef typename Superclass::EnergyDifferenceType EnergyDifferenceType;
    typedef typename Superclass::ExternalEnergyReturnType ExternalEnergyReturnType;
    typedef typename Superclass::DataImageType DataImageType;
    typedef typename Superclass::LabelImageType LabelImageType;
    typedef typename LabelImageType::OffsetType OffsetType;
    

//    typedef itk::Vector<EnergyDifferenceType, LabelImageType::ImageDimension> VectorType;
    typedef typename VectorImageType::PixelType VectorType;
    
    inline ExternalEnergyReturnType EvaluateEnergyDifference_(IndexType aInd, 
            LabelPixelType aLabelBefore, 
            LabelPixelType aLabelAfter,
            DataPixelType aImgValue);
   
    /** Add a point to the statistics of region. */
    inline void AddPoint_(IndexType aInd, LabelAbsPixelType aRegion, DataPixelType aVal);
    
    /** Remove a point from the statistics of a region. */
    inline void RemovePoint_(IndexType aInd /* = 0 */, LabelAbsPixelType aRegion, DataPixelType aVal);
    
    /** Implement the virtual method to clean-up the region statistics. */
    void KillRegion_(LabelAbsPixelType aL);
    
    /** Set the vector image to introduce a flow on the active contour. */
    void SetVectorImage(const VectorImageType *aPtr)
        { m_VecImage = aPtr; }
    
    /** Get the vector image that introduces a flow on the active contour*/
    const VectorImageType * GetDataImage() const
    { return m_VecImage.GetPointer(); }
    
    
    protected:
        RCVectorFieldFlow();
        ~RCVectorFieldFlow() {};
        void PrintSelf(std::ostream & os, Indent indent) const;
        

        
        
    private:
        RCVectorFieldFlow(const Self &);  //purposely not implemented
        void operator=(const Self &); //purposely not implemented
        
        /** Calculates a normalized vector in approximating the normal of the 
         * boundary of the object. The normal is returned via arg aRetvec and
         * the method returns the length of the normal.
         */
        EnergyDifferenceType ApproximateSurfaceNormalFromLabelImg(IndexType aIndex, 
                LabelAbsPixelType aL, VectorType* aRetvec);
        
                
        // energy related statistics statistics 
        typedef std::map<LabelAbsPixelType, EnergyDifferenceType> SumStatisticsMapType;
        typedef typename SumStatisticsMapType::iterator SumStatisticsIteratorType;
        
        /** Sum of intensities within a region - used to compute means */
        SumStatisticsMapType m_Sums;
        
        /** Sum of squares of intensities within regions - used to compute variances*/
        SumStatisticsMapType m_Sums_2;
        
        typedef std::map<LabelAbsPixelType, unsigned int>  CountStatisticsMapType;
        typedef typename CountStatisticsMapType::iterator CountStatisticsIteratorType;
        
        /** The number of pixel within a regin. */
        CountStatisticsMapType m_Count;
        
        VectorImageConstPointerType m_VecImage;
        

}; // end class
} // end namespace


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRCVectorFieldFlow.hxx"
#endif
 
#endif	/* ITKRCVECTORFIELDFLOW_H */

