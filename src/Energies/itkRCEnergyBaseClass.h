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
#ifndef __itkRCEnergyBaseClass_h
#define __itkRCEnergyBaseClass_h

#include "itkObject.h"
#include "itkObjectFactory.h"
//#include <type_traits>

namespace itk
{
/** \class 
 * \brief Base class for energy functions for the Region Competition optimizer and sampler.
 * 
 * \todo Should potentially be derived from itkFunctionBase.h or even itkImageFunction.h as
 * it is very similar to itkImageFunction. But the Evaluate methods do have a 
 * different interface.
 * 
 * \todo the LabelAbsPixelType should be an unsigned type of the LabelImage::PixelType.
 * std::make_unsigned in type_traits.h is experimental but could be used here.
 * Currently the type is set to unsigned int.
 * 
 * \ingroup Segmentation
 */
template <class TLabelImage, class TEnergyDifference = float>
class ITK_EXPORT RCEnergyBaseClass : public itk::Object
{
public:
  /** Standard class typedefs. */
  typedef RCEnergyBaseClass               Self;
  typedef Object                     Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** template related typedefs*/
  typedef TLabelImage LabelImageType;
  typedef typename LabelImageType::PixelType LabelPixelType;
  typedef typename LabelImageType::IndexType IndexType;
  typedef typename LabelImageType::RegionType RegionType;
  typedef typename LabelImageType::SizeType SizeType;
  typedef typename LabelImageType::OffsetType OffsetType;
  
  typedef unsigned int LabelAbsPixelType; // TODO!

  typedef TEnergyDifference EnergyDifferenceType;
    
  /** Run-time type information (and related methods). */
  itkTypeMacro(RCEnergyBaseClass, Object);

  /** Get the label image. */
  const LabelImageType * GetLabelImage() const
  { return m_LabelImage.GetPointer(); }
  
  /** Set the label image. */
  void SetLabelImage(const LabelImageType *aPtr)
  { m_LabelImage = aPtr; }
  
  /** Set the coefficient */
  itkSetMacro(Coefficient, EnergyDifferenceType);
  itkGetConstMacro(Coefficient, EnergyDifferenceType);

  /** Virtual function called by the optimizer after a region has been relabeled
   *  after a topological change. Instead of removing all the points this 
   *  is a shortcut to delete the region statistics.
   */
  virtual void KillRegion(LabelAbsPixelType aL) = 0;
  
  /** The function needs to be called once only before the evaluation is called 
   * first time. 
   * TODO: To avoid this energies should be part of the ITK pipeline.
   */
  virtual void PrepareEnergyCalculation(){};
  
  /**
   * The method is called before each iteration by the optimizer. It enables the
   * energy function to recompute statistics or images needed to compute 
   * the energy. For example a deconvolution energy might reconstruct a model 
   * image after each iteration or shape priors might update moment vectors.
   */
  virtual void PrepareEnergyCalculationForIteration(){};
  
  /** 
   * The method is called after having finished GenerateData(). This allows 
   * energies to clean up and release data structures. This is not the same as
   * destruction of the (energy) object. 
   */
  virtual void CleanUp(){};
  
  /** 
   * Calculate the total energy caused by this particular energy functional. 
   * Such calculations are usually very expensive and are used for debugging
   * or for plotting the energy over time. 
   */
  virtual EnergyDifferenceType CalculateTotalEnergy(){ return 0; };
  
protected:
  RCEnergyBaseClass();
  virtual ~RCEnergyBaseClass() {};
  void PrintSelf(std::ostream & os, Indent indent) const;

  // TODO: Outsource such helper methods to a utility class.
  
  /** Get the variance from the sum of squares, the mean, and the number of 
   * samples. A method frequently used by energies.*/
  inline EnergyDifferenceType CalculateVariance(EnergyDifferenceType aSumSq, 
          EnergyDifferenceType aMean, EnergyDifferenceType aN);
              
  /**
   * Get the sphere volume respecting the spacing of the label image.
   * A method commonly used by energies.
   */
  EnergyDifferenceType CalculateScaledSphereVolume(float aRadiusX);

  /** Const pointer to the input image. */
  typename LabelImageType::ConstPointer m_LabelImage;
  EnergyDifferenceType m_Coefficient;
  
private:
  RCEnergyBaseClass(const Self &);   //purposely not implemented
  void operator=(const Self &); //purposely not implemented


  
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRCEnergyBaseClass.hxx"
#endif

#endif
