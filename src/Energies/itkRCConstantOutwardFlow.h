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
#ifndef ITKRCCONSTANTOUTWARDFLOW_H
#define ITKRCCONSTANTOUTWARDFLOW_H

#include "itkRCInternalEnergyBaseClass.h"
#include <itkNumericTraits.h>


namespace itk
{

/**
 * @brief Returns a constant outward flow. See Cardinale et al. 2012.
 * 
 */
template<typename TLabelImage, typename TEnergyDifference = float>
class ITK_EXPORT RCConstantOutwardFlowEnergy :
  public RCInternalEnergyBaseClass<TLabelImage, TEnergyDifference>
  {
  public :
    /**
     * @name Standard ITK declarations
     */
    //@{
    typedef RCConstantOutwardFlowEnergy<TLabelImage, TEnergyDifference> Self;
    typedef RCInternalEnergyBaseClass<TLabelImage, TEnergyDifference> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<Self const> ConstPointer;

    itkNewMacro(Self);
    itkTypeMacro(RCConstantOutwardFlowEnergy, RCInternalEnergyBaseClass);
    
    typedef typename Superclass::IndexType IndexType;
    typedef typename Superclass::LabelPixelType LabelPixelType;
    typedef typename Superclass::LabelAbsPixelType LabelAbsPixelType;
    typedef typename Superclass::EnergyDifferenceType EnergyDifferenceType;
    typedef typename Superclass::InternalEnergyReturnType InternalEnergyReturnType;
    typedef typename Superclass::LabelImageType LabelImageType;
    typedef typename Superclass::OffsetType OffsetType;
    typedef typename Superclass::SizeType SizeType;
    
    inline InternalEnergyReturnType EvaluateEnergyDifference(IndexType aInd, 
            LabelPixelType aLabelBefore, 
            LabelPixelType aLabelAfter);
   
    /** Add a point to the statistics of region. */
    inline void AddPoint(IndexType aInd, LabelAbsPixelType aRegion);
    
    /** Remove a point from the statistics of a region. */
    inline void RemovePoint(IndexType aInd /* = 0 */, LabelAbsPixelType aRegion);
    
    /** Implement the virtual method to clean-up the region statistics. */
    void KillRegion(LabelAbsPixelType aL);
    
    protected:
        RCConstantOutwardFlowEnergy();
        ~RCConstantOutwardFlowEnergy() {};
        void PrintSelf(std::ostream & os, Indent indent) const;
        
    private:
        RCConstantOutwardFlowEnergy(const Self &);  //purposely not implemented
        void operator=(const Self &); //purposely not implemented
                    
        
}; // end class
} // end namespace


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRCConstantOutwardFlow.hxx"
#endif
 
#endif	/*  ITKRCCONSTANTOUTWARDFLOW_H */

