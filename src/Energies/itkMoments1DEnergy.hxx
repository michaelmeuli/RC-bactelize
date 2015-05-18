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
#ifndef _ITKMOMENTS1DDATACONTAINER_HXX
#define	_ITKMOMENTS1DDATACONTAINER_HXX

namespace itk{
    
    template <class TPoint, typename TStatistics>
    void Moments1DEnergy<TPoint, TStatistics>::
    Init(unsigned int aOrder, unsigned int aNbCentroids) {
        m_Order = aOrder;
        m_NbCentroids = (aNbCentroids < 1)? 1 : aNbCentroids;
        m_Coeff.resize(aNbCentroids);
        for(unsigned int vI = 0; vI < m_NbCentroids; vI++) {
            m_MomentsObjVec.push_back(Moments1DType::New());
            m_MomentsObjVec.back()->Init(m_Order, vI);
            m_MomentsObjVec.back()->SetPoints(&m_Points);
            m_Coeff[vI] = 1;
        }
        m_modified = true;
    }
    
    template <class TPoint, typename TStatistics>
    void Moments1DEnergy<TPoint, TStatistics>::
    AddPoint(PointType aPoint) {
        m_Points.insert(aPoint);
        m_modified = true;
    }
    
    template <class TPoint, typename TStatistics>
    void Moments1DEnergy<TPoint, TStatistics>::
    AddPoint(IndexType aPoint) {
        PointType vPoint;
        for(unsigned int vD = 0; vD < TPoint::GetPointDimension(); vD++){
            vPoint[vD] = aPoint[vD];
        }
        AddPoint(vPoint);
    }
    
    template <class TPoint, typename TStatistics>
    void Moments1DEnergy<TPoint, TStatistics>::
    RemovePoint(PointType aPoint) {
        m_Points.erase(aPoint);
        m_modified = true;
    }
    
    template <class TPoint, typename TStatistics>
    void Moments1DEnergy<TPoint, TStatistics>::
    RemovePoint(IndexType aPoint) {
        PointType vPoint;
        for(unsigned int vD = 0; vD < TPoint::GetPointDimension(); vD++){
            vPoint[vD] = aPoint[vD];
        }
        RemovePoint(vPoint);
    }
    
    template <class TPoint, typename TStatistics>
    void Moments1DEnergy<TPoint, TStatistics>::
    Update(){
        if(m_modified) {
            for(unsigned int vI = 0; vI < m_NbCentroids; vI++) {
                m_MomentsObjVec[vI]->Update();
            }
        }
        m_modified = false;
    }
 
    template <class TPoint, typename TStatistics>
    void Moments1DEnergy<TPoint, TStatistics>::
    SetCoefficient(unsigned int aIndex, float aCoefficient) {
        m_Coeff.at(aIndex) = aCoefficient;
    }
    
    template <class TPoint, typename TStatistics>
    typename Moments1DEnergy<TPoint, TStatistics>::StatisticsType 
    Moments1DEnergy<TPoint, TStatistics>::Dist1(Pointer aMomentsEnergy) {
        if(m_modified) {
            Update();
        }
        StatisticsType vDist1 = 0;
        for(unsigned int vI = 0; vI < m_NbCentroids; vI++) {
            vDist1 += m_Coeff[vI]*m_MomentsObjVec[vI]->Dist1(aMomentsEnergy->m_MomentsObjVec[vI]);
        }
        
    }
    
    template <class TPoint, typename TStatistics>
    typename Moments1DEnergy<TPoint, TStatistics>::StatisticsType 
    Moments1DEnergy<TPoint, TStatistics>::
    EnergyDiffL1Approx(const Pointer aMomentsEnergy, PointType aPoint, bool aAdd) {
        StatisticsType vEnergy = 0;
        typename Moments1DType::VectorType vUpdatedMoments;

        for(unsigned int vI = 0; vI < m_NbCentroids; vI++) {
            StatisticsType vDistBefore = aMomentsEnergy->m_MomentsObjVec[vI]
                    ->Dist1(this->m_MomentsObjVec[vI]);
            StatisticsType vDistAfter;
            if(aAdd) {
                m_MomentsObjVec[vI]->CalcSumOfPowersApproxWithPoint(aPoint, &vUpdatedMoments);
                vDistAfter = aMomentsEnergy->m_MomentsObjVec[vI]->Dist1(&vUpdatedMoments, m_Points.size() + 1);
            } else {
                m_MomentsObjVec[vI]->CalcSumOfPowersApproxWithoutPoint(aPoint, &vUpdatedMoments);
                vDistAfter = aMomentsEnergy->m_MomentsObjVec[vI]->Dist1(&vUpdatedMoments, m_Points.size() - 1);
            }
            vEnergy += m_Coeff[vI] * (vDistAfter - vDistBefore);
        }
        return vEnergy;
    }
}

#endif //_ITKMOMENTS1DDATACONTAINER_HXX
