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
#ifndef ITKRCDISCRETECONTOURLENGTHAPPROXENERGY_H
#define	ITKRCDISCRETECONTOURLENGTHAPPROXENERGY_H

#include "itkRCInternalEnergyBaseClass.h"
#include <itkNumericTraits.h>

#include "itkConstantBoundaryCondition.h"

namespace itk
{

/**
 * @brief Computes energy differences based on differences of local length 
 * measures using discrete approximations. 
 * 
 * For more details see "Y. Boykov and V. Kolmogorov, “Computing geodesics and minimal
 * surfaces via graph cuts,” in Proc. IEEE Int. Conf. Comput. Vis., vol. 1.
 * Nice, France, Oct. 2003, pp. 26–33.
 * 
 * A 16-neighborhood and 26-neighborhood relation is used to compute the 
 * geodesics as described in the paper above in 2D and 3D, resp.
 * 
 */
template<typename TLabelImage, typename TEnergyDifference = float>
class ITK_EXPORT RCDiscreteContourLengthApproxEnergyEnergy :
  public RCInternalEnergyBaseClass<TLabelImage, TEnergyDifference>
  {
  public :
    /**
     * @name Standard ITK declarations
     */
    //@{
    typedef RCDiscreteContourLengthApproxEnergyEnergy<TLabelImage, TEnergyDifference> Self;
    typedef RCInternalEnergyBaseClass<TLabelImage, TEnergyDifference> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<Self const> ConstPointer;

    itkNewMacro(Self);
    itkTypeMacro(RCDiscreteContourLengthApproxEnergyEnergy, RCInternalEnergyBaseClass);
    
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
        RCDiscreteContourLengthApproxEnergyEnergy();
        ~RCDiscreteContourLengthApproxEnergyEnergy() {};
        void PrintSelf(std::ostream & os, Indent indent) const;
        
        InternalEnergyReturnType EvaluateAtIndexUsingLabel(IndexType const & aIndex, LabelPixelType aL) const;
        
    private:
        RCDiscreteContourLengthApproxEnergyEnergy(const Self &);  //purposely not implemented
        void operator=(const Self &); //purposely not implemented
                    
        static inline int H(LabelPixelType aX){ return (aX > 0) ? 1 : 0; }
        
        unsigned int m_Radius;
        std::vector<float> m_lengthWeightMask;
        std::vector<OffsetType> m_Offsets;

        bool m_NeedToUseBoundaryCondition;    
        
}; // end class
} // end namespace


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRCDiscreteContourLengthApproxEnergy.hxx"
#endif
 
#endif	/* ITKRCDISCRETECONTOURLENGTHAPPROXENERGY_H */

