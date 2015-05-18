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
#ifndef ITKRCPIECEWISESMOOTHLIENERGY_H
#define	ITKRCPIECEWISESMOOTHLIENERGY_H

#include "itkRCExternalEnergyBaseClass.h"
#include "itkSphereBitmapImageSource.h"
#include "itkGaussianImageSource.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include <map>


namespace itk
{

/**
 * @brief Computes energy differences for a pixel change in for a piece-wise 
 * smooth image model and with i.i.d. Gaussian noise. 
 *
 * The energy to minimize is ...
 * 
 * 
 */
template<typename TLabelImage, typename TDataImage, typename TEnergyDifference = float>
class ITK_EXPORT RCPiecewiseSmoothLiEnergy :
  public RCExternalEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference>
  {
  public :
    /**
     * @name Standard ITK declarations
     */
    //@{
    typedef RCPiecewiseSmoothLiEnergy<TLabelImage, TDataImage, TEnergyDifference> Self;
    typedef RCExternalEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<Self const> ConstPointer;

    itkNewMacro(Self);
    itkTypeMacro(RCPiecewiseSmoothLiEnergy, RCExternalEnergyBaseClass);
    
    typedef typename Superclass::IndexType IndexType;
    typedef typename Superclass::LabelPixelType LabelPixelType;
    typedef typename Superclass::LabelAbsPixelType LabelAbsPixelType;
    typedef typename Superclass::LabelImageType LabelImageType;
    typedef typename Superclass::DataPixelType DataPixelType;
    typedef typename Superclass::DataImageType DataImageType;
    typedef typename Superclass::EnergyDifferenceType EnergyDifferenceType;
    typedef typename Superclass::ExternalEnergyReturnType ExternalEnergyReturnType;


    
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
    
    /** Implements the virtual method. */
    void PrepareEnergyCalculation();
    
    itkSetMacro(Sigma, float);
    itkGetConstMacro(Sigma, float);
    
    protected:
        RCPiecewiseSmoothLiEnergy();
        ~RCPiecewiseSmoothLiEnergy() {};
        void PrintSelf(std::ostream & os, Indent indent) const;
        
    private:
        RCPiecewiseSmoothLiEnergy(const Self &);  //purposely not implemented
        void operator=(const Self &); //purposely not implemented
        
        // 
        typedef GaussianImageSource<DataImageType> GaussianImageSourceType;
        typename GaussianImageSourceType::Pointer m_GaussianImageSource;
        
        
        // energy related members
        typedef bool MaskPixelType;
        typedef Image<MaskPixelType, LabelImageType::ImageDimension> MaskImageType;
        typedef SphereBitmapImageSource<MaskImageType> SphereMaskImageSourceType;
        typedef typename SphereMaskImageSourceType::Pointer SphereMaskImagePointerType;
        SphereMaskImagePointerType m_SphereMaskForLocalEnergy;
        
        float m_Sigma;
                

        
}; // end class
} // end namespace


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRCPiecewiseSmoothLiEnergy.hxx"
#endif
 
#endif	/* ITKRCPIECEWISESMOOTHLIENERGY_H */

