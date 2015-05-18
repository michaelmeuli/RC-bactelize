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
#ifndef ITKRCGPUENERGYBASECLASS_H
#define	ITKRCGPUENERGYBASECLASS_H

#include <vector>
#include <utility>
#include "itkRCEnergyBaseClass.h"


namespace itk {

    /** \class 
     * \brief Abstract base class for energy functions using the GPU. Compared 
     * to internal and external energies, this classes derived from here 
     * incorporate both, the external and the internal term. Hence, usually only
     * one energy class is registered at the (RC) optimizer.
     * Here, all energy differences are computed in one shot. Therefor there are 
     * no AddPoint or RemovePoint methods. 
     * 
     * \ingroup Segmentation
     */
    template <class TLabelImage, class TDataImage, class TEnergyDifference = float>
    class ITK_EXPORT RCGPUEnergyBaseClass
    : public RCEnergyBaseClass<TLabelImage, TEnergyDifference> {
    public:
        /** Standard class typedefs. */
        typedef RCGPUEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference> Self;
        typedef RCEnergyBaseClass<TLabelImage, TEnergyDifference> Superclass;
        typedef SmartPointer< Self > Pointer;
        typedef SmartPointer< const Self > ConstPointer;
        
        itkTypeMacro(RCGPUEnergyBaseClass, RCEnergyBaseClass);

        /** derived typedefs*/
        typedef typename Superclass::LabelPixelType LabelPixelType;
        typedef typename Superclass::LabelAbsPixelType LabelAbsPixelType;
        typedef typename Superclass::LabelImageType LabelImageType;
        typedef typename Superclass::IndexType IndexType;
        typedef typename Superclass::RegionType RegionType;
        typedef typename Superclass::SizeType SizeType; 
        typedef typename Superclass::OffsetType OffsetType;
                                
        typedef typename Superclass::EnergyDifferenceType EnergyDifferenceType;        

        /** Data image typedefs*/
        typedef TDataImage DataImageType;
        typedef typename DataImageType::PixelType DataPixelType;
        
        /* Begin of GPU typedefs */
        typedef int32_t GPU_IntType;
        typedef uint32_t GPU_UIntType;

          /** Get the data image. */
        const DataImageType * GetDataImage() const
        { return m_DataImage.GetPointer(); }
  
        /** Set the label image. This needs to be called to define a valid 
         * external energy.
         */
        void SetDataImage(const DataImageType *aPtr)
        { m_DataImage = aPtr; }
  
        /** Add a point to the statistics of region. */
        virtual void AddPoint(IndexType aInd, LabelAbsPixelType aLabel, DataPixelType aVal) = 0;
        
        /** Remove a point from the statistics of a region. */
        virtual void RemovePoint(IndexType aInd /* = 0 */, LabelAbsPixelType aLabel, DataPixelType aVal) = 0;
        
        
        
        /** The default behaviour is to call AddPoint(...) with the label where the
         *  point will move to and RemovePoint(...) with the label of which the point
         *  is "coming" from. For some energies it might be beneficial to 
         *  overwrite this method.
         */
        virtual void SwitchPoint(IndexType aInd, LabelAbsPixelType aLabelFrom, 
                LabelAbsPixelType aLabelTo, DataPixelType aVal) {
            AddPoint(aInd, aLabelTo, aVal);
            RemovePoint(aInd, aLabelFrom, aVal);
        };
        
  
        /**
         * Batch processing for energy difference evaluation. The method takes
         * references to vectors in order to avoid frequent function calls to 
         * virtual method EvaluateEnergyDifference. 
         */
        virtual void EvaluateEnergyDifferences(unsigned int **aGPUCoordinates,
                unsigned int* aGPUToLabel,
                float* aGPUOutInternalEnergy,
                float* aGPUOutExternalEnergy,
                unsigned int* aGPUOutMerge,
                unsigned int aPaddedNbCandidates
        ) = 0;

        virtual void PrepareEnergyCalculation(){};
        virtual void PrepareEnergyCalculationForIteration(){};
        
        protected:
            RCGPUEnergyBaseClass();
            ~RCGPUEnergyBaseClass() {};
            void PrintSelf(std::ostream & os, Indent indent) const;
            
            /** The data image is needed for external energies. */
            typedef typename DataImageType::ConstPointer DataImageConstPointerType;
            DataImageConstPointerType m_DataImage;
            

        private:
            RCGPUEnergyBaseClass(const Self &);  //purposely not implemented
            void operator=(const Self &); //purposely not implemented
            
        public:

            
    }; // class

} // namespace


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRCGPUEnergyBaseClass.hxx"
#endif

#endif	/* ITKRCGPUENERGYBASECLASS_H */

