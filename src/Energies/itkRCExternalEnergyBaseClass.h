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
#ifndef ITKRCEXTERNALENERGYBASECLASS_H
#define	ITKRCEXTERNALENERGYBASECLASS_H

#include <vector>
#include <utility>
#include "itkRCEnergyBaseClass.h"


namespace itk {

    /** \class 
     * \brief Base class for energy functions for the Region Competition optimizer and sampler.
     *
     * \ingroup Segmentation
     */
    template <class TLabelImage, class TDataImage, class TEnergyDifference = float>
    class ITK_EXPORT RCExternalEnergyBaseClass
    : public RCEnergyBaseClass<TLabelImage, TEnergyDifference> {
    public:
        /** Standard class typedefs. */
        typedef RCExternalEnergyBaseClass<TLabelImage, TDataImage, TEnergyDifference> Self;
        typedef RCEnergyBaseClass<TLabelImage, TEnergyDifference> Superclass;
        typedef SmartPointer< Self > Pointer;
        typedef SmartPointer< const Self > ConstPointer;
        
        itkTypeMacro(RCExternalEnergyBaseClass, RCEnergyBaseClass);

        /** derived typedefs*/
        typedef typename Superclass::LabelPixelType LabelPixelType;
        typedef typename Superclass::LabelAbsPixelType LabelAbsPixelType;
        typedef typename Superclass::IndexType IndexType;
        typedef typename Superclass::RegionType RegionType;
        typedef typename Superclass::SizeType SizeType; 
        typedef typename Superclass::OffsetType OffsetType;
                                
        typedef typename Superclass::EnergyDifferenceType EnergyDifferenceType;        

        /** Data image typedefs*/
        typedef TDataImage DataImageType;
        typedef typename DataImageType::PixelType DataPixelType;

        /** Specific typedefs*/
        typedef std::pair<EnergyDifferenceType, bool> ExternalEnergyReturnType;

          /** Get the data image. */
        const DataImageType * GetDataImage() const
        { return m_DataImage.GetPointer(); }
  
        /** Set the label image. This needs to be called to define a valid 
         * external energy.
         */
        void SetDataImage(const DataImageType *aPtr)
        { m_DataImage = aPtr; }

        public:
        /** Add a point to the statistics of region. */
        virtual void AddPoint(IndexType aInd, LabelAbsPixelType aLabel, DataPixelType aVal);
        
        protected:
        /** Add a point to the statistics of a region - derived class implementation*/
        virtual void AddPoint_(IndexType aInd, LabelAbsPixelType aLabel, DataPixelType aVal) = 0;
        
        public:
        /** Remove a point from the statistics of a region. */
        virtual void RemovePoint(IndexType aInd /* = 0 */, LabelAbsPixelType aLabel, DataPixelType aVal);

        protected:
        /** Remove a point to the statistics of a region - derived class implementation*/
        virtual void RemovePoint_(IndexType aInd /* = 0 */, LabelAbsPixelType aLabel, DataPixelType aVal) = 0;
        
        public:
        virtual void KillRegion(LabelAbsPixelType aL);
            
        protected:
        virtual void KillRegion_(LabelAbsPixelType aL) = 0;
                
        public:
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
        
        /** This method returns the difference in energy when changing the 
         *   pixel with index aInd from aLabelBefore to aLabelAfter.
         */   
        virtual ExternalEnergyReturnType EvaluateEnergyDifference(IndexType aInd, 
                LabelPixelType aLabelBefore, 
                LabelPixelType aLabelAfter,
                DataPixelType aImgValue);
        
        protected:
        virtual ExternalEnergyReturnType EvaluateEnergyDifference_(IndexType aInd, 
                LabelPixelType aLabelBefore, 
                LabelPixelType aLabelAfter,
                DataPixelType aImgValue) = 0;
        
        public:        
        /**
         * Batch processing for energy difference evaluation. The method takes
         * references to vectors in order to avoid frequent function calls to 
         * virtual method EvaluateEnergyDifference. The last argument contains 
         * the result.
         * @see EvaluateEnergyDifference
         */
        virtual void EvaluateEnergyDifferences(const std::vector<IndexType>& aInds, 
                const std::vector<LabelPixelType>& aLabelBefore, 
                const std::vector<LabelPixelType>& aLabelAfter,
                const std::vector<DataPixelType>& aImgValue,
                std::vector<ExternalEnergyReturnType>& aReturnVec);
        
        typedef DataImageType ReconstructedImageType;
        /** 
         * The method generates (if implemented) a reconstructed image based on 
         * the current label image and the energies image formation model.
         * The default behaviour throws an exception.
         * 
         * Usage: 
         * typedef EnergyType::ReconstructedImageType ReconstructedImageType;
         * ReconstructedImageType::Pointer vRecImg = ReconstructedImageType::New();
         * void * vimg = &vRecImg;     // cast to a void pointer
         * GenerateModelImage(vimg);   // pass the void pointer
         * 
         */
        virtual void GenerateReconstructedImage(void* aPointerToResultImage) {
            itkExceptionMacro("The method is not implemented for this energy. ");
        };
        
        itkSetMacro(RegionMergingThreshold, EnergyDifferenceType);
        itkGetConstMacro(RegionMergingThreshold, EnergyDifferenceType);

        protected:
            RCExternalEnergyBaseClass();
            ~RCExternalEnergyBaseClass() {};
            void PrintSelf(std::ostream & os, Indent indent) const;
            
            /** The data image is needed for external energies. */
            typedef typename DataImageType::ConstPointer DataImageConstPointerType;
            DataImageConstPointerType m_DataImage;
            EnergyDifferenceType m_RegionMergingThreshold;
            

        private:
            RCExternalEnergyBaseClass(const Self &);  //purposely not implemented
            void operator=(const Self &); //purposely not implemented
            
        public:
            // todo: should be static
            inline EnergyDifferenceType CalculateKullbackLeiblerDistance(
            EnergyDifferenceType aMu1, EnergyDifferenceType aMu2, 
                    EnergyDifferenceType aVar1,EnergyDifferenceType aVar2,
                    EnergyDifferenceType aN1, EnergyDifferenceType aN2);
            
        protected:
            
            /** Standard statistics (used by almost any external energy). */
            
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
            
    }; // class

} // namespace


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRCExternalEnergyBaseClass.hxx"
#endif

#endif	/* ITKRCEXTERNALENERGYBASECLASS_H */

