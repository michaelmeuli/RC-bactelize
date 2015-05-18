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

#include "itkRCShapePrior1DEnergy.h"
#include "itkSphereBitmapImageSource.h"
#include "itkRCEnergyBaseClass.h"



namespace itk
{
    template<typename TLabelImage, typename TTemplateImage, typename TEnergyDifference>
    RCShapePrior1DEnergy<TLabelImage, TTemplateImage, TEnergyDifference>
    ::RCShapePrior1DEnergy() {
        m_BackgroundValue = 0;
        m_BackgroundLabel = 0;
        m_OrderOfMoments = 20;
        m_ShapePriorMinSize = 50;
        m_ShapePriorNbCentroids = 1;
        // the coefficients will be reset in Init()/PrepareEnergyCalculation();
        m_ShapePriorCentroidCoeff.Fill(1);
        m_Moments1DReference = Moments1DType::New();
    }
    
    template<typename TLabelImage, typename TTemplateImage, typename TEnergyDifference>
    void RCShapePrior1DEnergy<TLabelImage, TTemplateImage, TEnergyDifference>
    ::PrintSelf(std::ostream & os, Indent indent) const
    {
        Superclass::PrintSelf(os, indent);
        os << indent << "Background value: " << m_BackgroundValue <<  std::endl;
        os << indent << "Background label:" << m_BackgroundLabel << std::endl;
        os << indent << "Order of moments: " << m_OrderOfMoments << std::endl; 
        os << indent << "Minimum size to bias with shape prior: " << m_ShapePriorMinSize << std::endl;
        os << indent << "Number of centroids and their coefficients: " << m_ShapePriorNbCentroids << std::endl;
        for(unsigned int vI = 0; vI < m_ShapePriorNbCentroids; vI++) {
            os << indent << "Coefficient [" << vI << "]: " << m_ShapePriorCentroidCoeff << std::endl;
        }
    }
    
    template<typename TLabelImage, typename TTemplateImage, typename TEnergyDifference>
    void RCShapePrior1DEnergy<TLabelImage, TTemplateImage, TEnergyDifference>
    ::PrepareEnergyCalculationForIteration() {
        
            Moments1DStatisticsIteratorType vActiveLabelsIt =
                    m_Moments1DStatistics.begin();
            Moments1DStatisticsIteratorType vActiveLabelsItEnd =
                    m_Moments1DStatistics.end();

            for (; vActiveLabelsIt != vActiveLabelsItEnd; ++vActiveLabelsIt) {
                LabelAbsPixelType vLabelAbs = vActiveLabelsIt->first;
                if (vLabelAbs == 0) {
                    continue;
                }
                m_Moments1DStatistics[vLabelAbs]->Update();
            }
    }
    template<typename TLabelImage, typename TTemplateImage, typename TEnergyDifference>
    void RCShapePrior1DEnergy<TLabelImage, TTemplateImage, TEnergyDifference>
    ::PrepareEnergyCalculation() {
        
//        typename Moments1DType::Pointer vMomentsMutant = Moments1DType::New();

        itk::ImageRegionConstIteratorWithIndex<TemplateImageType> vIt(m_TemplateImage, 
                m_TemplateImage->GetBufferedRegion());

        for (vIt.GoToBegin(); !vIt.IsAtEnd(); ++vIt) {
            if (vIt.Get() == m_BackgroundValue) {
                continue;
            }

            // convert the index to a point:
            typename Moments1DType::PointType vPoint;
            typename TemplateImageType::IndexType vIndex = vIt.GetIndex();
            for (int vD = 0; vD < TemplateImageType::ImageDimension; vD++) {
                vPoint[vD] = vIndex[vD];
            }
            m_Moments1DReference->AddPoint(vPoint);
//            vMomentsMutant->AddPoint(vPoint);
        }

        m_Moments1DReference->Init(m_OrderOfMoments, m_ShapePriorNbCentroids);
        //m_Moments1DReference->SetNbCentroids(m_ShapePriorNbCentroids);
        for(unsigned int vI = 0; vI < m_ShapePriorNbCentroids; vI++) {
            m_Moments1DReference->SetCoefficient(vI,m_ShapePriorCentroidCoeff[vI]);
        }
        
        m_Moments1DReference->Update();
    
    }
        
    template<typename TLabelImage, typename TTemplateImage, typename TEnergyDifference>
    inline 
    typename RCShapePrior1DEnergy<TLabelImage, TTemplateImage, TEnergyDifference>::
    InternalEnergyReturnType
    RCShapePrior1DEnergy<TLabelImage, TTemplateImage, TEnergyDifference>
    ::EvaluateEnergyDifference(IndexType aInd, 
            LabelPixelType aLabelFrom, 
            LabelPixelType aLabelTo) {
        
        InternalEnergyReturnType vDeltaEnergy = 0;
        typename Moments1DType::PointType vPoint;
        for (unsigned int vD = 0; vD < TemplateImageType::ImageDimension; vD++) {
            vPoint[vD] = aInd[vD];
        }

        if (aLabelFrom != 0 && 
                m_Moments1DStatistics[aLabelFrom]->GetSize() >= m_ShapePriorMinSize) { // shrinking or competition
            
            vDeltaEnergy += m_Moments1DStatistics[aLabelFrom]->EnergyDiffL1Approx(m_Moments1DReference, vPoint, false);
        }
        
        if (aLabelTo != 0 && 
                m_Moments1DStatistics[aLabelTo]->GetSize() >= m_ShapePriorMinSize) { // growing or competition
            
            vDeltaEnergy += m_Moments1DStatistics[aLabelTo]->EnergyDiffL1Approx(m_Moments1DReference, vPoint, true);
        }

        return this->m_Coefficient * vDeltaEnergy;
    }
    
    template<typename TLabelImage, typename TTemplateImage, typename TEnergyDifference>
    inline void 
    RCShapePrior1DEnergy<TLabelImage, TTemplateImage, TEnergyDifference>
    ::AddPoint(IndexType aInd, LabelAbsPixelType aLabel) {
        
        if(aLabel == m_BackgroundLabel) return;
        
        Moments1DStatisticsIteratorType vIt = m_Moments1DStatistics.find(aLabel);
        if(vIt == m_Moments1DStatistics.end()) {
            m_Moments1DStatistics[aLabel] = Moments1DType::New();
            m_Moments1DStatistics[aLabel]->Init(m_OrderOfMoments, m_ShapePriorNbCentroids);
            for(unsigned int vI = 0; vI < m_ShapePriorNbCentroids; vI++) {
                m_Moments1DStatistics[aLabel]->SetCoefficient(vI,m_ShapePriorCentroidCoeff[vI]);
            }
        }
        
        m_Moments1DStatistics[aLabel]->AddPoint(aInd);
    }
    
    template<typename TLabelImage, typename TTemplateImage, typename TEnergyDifference>
    inline void 
    RCShapePrior1DEnergy<TLabelImage, TTemplateImage, TEnergyDifference>
    ::RemovePoint(IndexType aInd, LabelAbsPixelType aLabel) {
        
        if(aLabel == m_BackgroundLabel) return;
        
        assert(m_Moments1DStatistics.find(aLabel) != m_Moments1DStatistics.end());
        if(m_Moments1DStatistics[aLabel]->GetSize() == 1) {
            KillRegion(aLabel);
        } else {
            m_Moments1DStatistics[aLabel]->RemovePoint(aInd);
        }
    }
    
    template<typename TLabelImage, typename TTemplateImage, typename TEnergyDifference>
    void 
    RCShapePrior1DEnergy<TLabelImage, TTemplateImage, TEnergyDifference>
    ::KillRegion(LabelAbsPixelType aLabel) {
        m_Moments1DStatistics.erase(aLabel);
    }
        
} // end namespace itk