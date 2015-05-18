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
#ifndef ITKRCPIECEWISECONSTANTSQUAREDISTNERGY_H
#define	ITKRCPIECEWISECONSTANTSQUAREDISTNERGY_H

#include "itkRCExternalEnergyBaseClass.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include <map>
namespace itk
{

/**
 * @brief Computes energy differences for a pixel change in for a piece-wise 
 * constant image model and with i.i.d. Gaussian noise. 
 *
 * The energy to minimize is \f$E = \sum_i^M(\mu_i - I(x))^2\f$ with M being the 
 * number of regions, I the image and \f$\mu_i\f$ the mean of region i.
 * 
 * In case of 2 regions this image model is also called Chan-Vese model. 
 * 
 */
template<typename TLabelImage, typename TDataImage, typename TEnergyDifference = float>
class ITK_EXPORT RCPiecewiseConstantSquareDistEnergy :
  public RCExternalEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
  {
  public :
    /**
     * @name Standard ITK declarations
     */
    //@{
    typedef RCPiecewiseConstantSquareDistEnergy<TLabelImage, TDataImage, TEnergyDifference> Self;
    typedef RCExternalEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<Self const> ConstPointer;

    itkNewMacro(Self);
    itkTypeMacro(RCPiecewiseConstantSquareDistEnergy, RCExternalEnergyBaseClass);
    
    typedef typename Superclass::IndexType IndexType;
    typedef typename Superclass::LabelPixelType LabelPixelType;
    typedef typename Superclass::LabelImageType LabelImageType;
    typedef typename Superclass::LabelAbsPixelType LabelAbsPixelType;
    typedef typename Superclass::DataImageType DataImageType;
    typedef typename Superclass::DataPixelType DataPixelType;
    typedef typename Superclass::EnergyDifferenceType EnergyDifferenceType;
    typedef typename Superclass::ExternalEnergyReturnType ExternalEnergyReturnType;
    typedef typename Superclass::ReconstructedImageType ReconstructedImageType;

    
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
    
    /** Overrides the base class method */
    EnergyDifferenceType CalculateTotalEnergy();
    
    /** Overrides the base class method */
    void GenerateReconstructedImage(void* aPointerToResultImage);
    
    protected:
        RCPiecewiseConstantSquareDistEnergy();
        ~RCPiecewiseConstantSquareDistEnergy() {};
        void PrintSelf(std::ostream & os, Indent indent) const;
        
        
    private:
        RCPiecewiseConstantSquareDistEnergy(const Self &);  //purposely not implemented
        void operator=(const Self &); //purposely not implemented
                

        
}; // end class
} // end namespace


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRCPiecewiseConstantSquareDistEnergy.hxx"
#endif
 
#endif	/* ITKRCPIECEWISECONSTANTSQUAREDISTNERGY_H */

