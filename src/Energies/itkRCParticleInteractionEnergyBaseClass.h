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
#ifndef ITKRCPARTICLEINTERACTIONENERGYBASECLASS_H
#define	ITKRCPARTICLEINTERACTIONENERGYBASECLASS_H

#include <vector>
#include <utility>
#include "itkRCEnergyBaseClass.h"
#include "itkConnectivity.h"
#include "itkBackgroundConnectivity.h"
#include "itkVector.h"

namespace itk {

    /** \class 
     * \brief Baseclass for (internal) energies that involve particle interactions 
     * for the Region competition framework. The speciality is that the derived
     * classes have access to the particle container of the optimizer. 
     *
     * \ingroup Segmentation
     */
    template <class TLabelImage, class TOptimizer, class TEnergyDifference = float>
    class ITK_EXPORT RCParticleInteractionEnergyBaseClass
    : public RCEnergyBaseClass<TLabelImage, TEnergyDifference> {
    public:
        /** Standard class typedefs. */
        typedef RCParticleInteractionEnergyBaseClass<TLabelImage, TOptimizer, TEnergyDifference> Self;
        typedef RCEnergyBaseClass<TLabelImage, TEnergyDifference> Superclass;
        typedef SmartPointer< Self > Pointer;
        typedef SmartPointer< const Self > ConstPointer;

        typedef TOptimizer OptimizerType;
        typedef typename OptimizerType::InteractionParticlesIteratorType InteractionParticlesIteratorType;
        
        itkTypeMacro(RCParticleInteractionEnergyBaseClass, RCEnergyBaseClass);

        /** derived typedefs*/
        typedef typename Superclass::LabelPixelType LabelPixelType;
        typedef typename Superclass::LabelAbsPixelType LabelAbsPixelType;
        typedef typename Superclass::LabelImageType LabelImageType;
        typedef typename Superclass::IndexType IndexType;
        typedef typename Superclass::RegionType RegionType;
        typedef typename Superclass::SizeType SizeType; 
        typedef typename Superclass::OffsetType OffsetType;
                                
        typedef typename Superclass::EnergyDifferenceType EnergyDifferenceType;        


        /** Specific typedefs*/
        typedef EnergyDifferenceType ParticleInteractionEnergyReturnType;
  
        /** Add a point to the statistics of region. */
        virtual void AddPoint(IndexType aInd, LabelAbsPixelType aLabel){};
        
        /** Remove a point from the statistics of a region. */
        virtual void RemovePoint(IndexType aInd, LabelAbsPixelType aLabel){};
  
        /** The default behaviour is to call AddPoint(...) with the label where the
         *  point will move to and RemovePoint(...) with the label of which the point
         *  is "coming" from. For some energies it might be beneficial to 
         *  overwrite this method.
         */
        virtual void SwitchPoint(IndexType aInd, LabelAbsPixelType aLabelFrom, 
                LabelAbsPixelType aLabelTo) {
            AddPoint(aInd, aLabelTo);
            RemovePoint(aInd, aLabelFrom);
        };
        
        /** This method returns the regularizing flow for a the current label.
         */   
        virtual ParticleInteractionEnergyReturnType EvaluateEnergyDifference(IndexType aInd, 
                InteractionParticlesIteratorType aBegin, InteractionParticlesIteratorType aEnd) = 0;
        
        /**
         * Batch processing for energy difference evaluation. The method takes
         * references to vectors in order to avoid frequent function calls to 
         * virtual method EvaluateEnergyDifference. The last argument contains 
         * the result.
         * @see EvaluateEnergyDifference
         */
//        virtual void EvaluateEnergyDifferences(const std::vector<IndexType>& aInds, 
//                const std::vector<LabelPixelType>& aLabelBefore, 
//                const std::vector<LabelPixelType>& aLabelAfter,
//                std::vector<ParticleInteractionEnergyReturnType>& aReturnVec);
        
            
        itkSetMacro(ParticleInteractionRadius, float);
        itkGetMacro(ParticleInteractionRadius, float);
        
        virtual SizeType GetMinimalIteratorSize() = 0;
        
        protected:
            RCParticleInteractionEnergyBaseClass();
            ~RCParticleInteractionEnergyBaseClass() {};
            void PrintSelf(std::ostream & os, Indent indent) const;
            
            typedef itk::Vector<EnergyDifferenceType, LabelImageType::ImageDimension> VectorType;
            /** Calculates a normalized vector in approximating the normal of the 
             * boundary of the object. The normal is returned via arg aRetvec and
             * the method returns the length of the normal.
             */
            virtual EnergyDifferenceType ApproximateSurfaceNormal(IndexType aIndex, VectorType& aRetvec);
            
            
        private:
            RCParticleInteractionEnergyBaseClass(const Self &);  //purposely not implemented
            void operator=(const Self &); //purposely not implemented

        protected:
            typedef Connectivity<LabelImageType::ImageDimension, LabelImageType::ImageDimension - 1 > 
            ForegroundConnectivityType;
            
            typedef typename BackgroundConnectivity<ForegroundConnectivityType>::Type 
            BackgroundConnectivityType;
        
            unsigned int m_NeighborhoodSize_BG_Connectivity;
            const OffsetType * m_NeighborsOffsets_BG_Connectivity;
            
            float m_ParticleInteractionRadius;
            

        public:
            
    }; // class

} // namespace


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRCParticleInteractionEnergyBaseClass.hxx"
#endif

#endif	/* ITKRCPARTICLEINTERACTIONENERGYBASECLASS_H */

