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
#ifndef _ITKMOMENTS1D_HXX
#define	_ITKMOMENTS1D_HXX

#include <set>
#include <utility>
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "itkMoments1D.h"

namespace itk {

    template <class TPoint, class TStatistics> 
    bool Moments1D<TPoint, TStatistics>::m_NchooseK_initialized = false;
    
    template <class TPoint, class TStatistics> 
    typename Moments1D<TPoint, TStatistics>::MatrixType 
    Moments1D<TPoint, TStatistics>::m_NChooseK;

    template <class TPoint, class TStatistics>
    void Moments1D<TPoint, TStatistics>::
    Init(unsigned int aOrder) {
        Init(aOrder,0);
    }
     
    template <class TPoint, class TStatistics>
    void Moments1D<TPoint, TStatistics>::
    Init(unsigned int aOrder, unsigned char aType) {
        m_CentroidType = aType;
        m_SumOfPowers.set_size(aOrder + 1);
        if(!Self::m_NchooseK_initialized || aOrder > m_Order) {
            initializeNchooseK(aOrder);
            Self::m_NchooseK_initialized = true;
        }
        m_Order = aOrder;
    }

    template <class TPoint, class TStatistics>
    unsigned int Moments1D<TPoint, TStatistics>::
    GetSize() const {
        return m_Points->size();
    }

    template <class TPoint, class TStatistics>
    void Moments1D<TPoint, TStatistics>::
    CalcSumOfPowersWithPoint(TPoint aPoint, VectorType* aSumOfPowersOutput)  {

        std::pair<PointContainerIteratorType, bool> vRet;
        vRet = m_Points->insert(aPoint);

        PointType vCOM = CalcCenterOfMass(m_Points);
        
        PointType vCentroid = CalcCentroid(m_Points, vCOM);

        StatisticsType vMu = CalcMu(m_Points, vCentroid);

        aSumOfPowersOutput->set_size(m_Order+1);
        CalcSumOfPowers(m_Points, vCentroid, vMu, aSumOfPowersOutput);

        m_Points->erase(vRet.first);
    }

    template <class TPoint, class TStatistics>
    void Moments1D<TPoint, TStatistics>::
    CalcSumOfPowersWithoutPoint(TPoint aPoint, VectorType* aSumOfPowersOutput)  {

        m_Points->erase(aPoint);

        PointType vCOM = CalcCenterOfMass(m_Points);

        PointType vCentroid = CalcCentroid(m_Points, vCOM);
                
        StatisticsType vMu = CalcMu(m_Points, vCentroid);

        aSumOfPowersOutput->set_size(m_Order+1);
        CalcSumOfPowers(m_Points, vCentroid, vMu, aSumOfPowersOutput);

        m_Points->insert(aPoint);
    }

    /// This method approximates the updated moments when adding a multidimensional
    /// point to the data set. It is an approximation because the 1D moments
    /// are based on distances referring to the center of mass of the point set 
    /// but the center of mass here is not updated.
    template <class TPoint, class TStatistics>
    void Moments1D<TPoint, TStatistics>::
    CalcSumOfPowersApproxWithPoint(TPoint aPoint, VectorType* aSumOfPowersOutput) const{
        
        aSumOfPowersOutput->set_size(m_Order+1);

        double vNewDist = aPoint.EuclideanDistanceTo(m_Centroid);

        UpdateSumOfPowersWithData(&m_SumOfPowers, aSumOfPowersOutput, m_Points->size(), m_Mu, vNewDist);
    }

 
    /// This method approximates the updated moments when adding a multidimensional
    /// point to the data set. It is an approximation because the 1D moments
    /// are based on distances referring to the center of mass of the point set 
    /// but the center of mass here is not updated.
    template <class TPoint, class TStatistics>
    void Moments1D<TPoint, TStatistics>::
    CalcSumOfPowersApproxWithoutPoint(TPoint aPoint, VectorType* aSumOfPowersOutput) const{

        unsigned int vN = this->GetSize();

        // update the center of mass first (to have a symmetric metric)
//        PointType vCOM;
//        for(unsigned int vD = 0; vD < PointType::Dimension; vD++) {
//            vCOM[vD] = (m_CenterOfMass[vD] * vN - aPoint[vD]) / (vN - 1);
//        }

        double vOldDist = aPoint.EuclideanDistanceTo(m_Centroid);

        aSumOfPowersOutput->set_size(m_Order + 1);
        UpdateSumOfPowersWithoutData(&m_SumOfPowers, aSumOfPowersOutput, vN, m_Mu, vOldDist);
    }


    template <class TPoint, class TStatistics>
    typename Moments1D<TPoint, TStatistics>::StatisticsType
    Moments1D<TPoint, TStatistics>::
    Dist1(const Pointer aM) const {
        return Dist1(&(this->m_SumOfPowers), this->GetSize(),
                &(aM->m_SumOfPowers), aM->GetSize());
    }

    template <class TPoint, class TStatistics>
    typename Moments1D<TPoint, TStatistics>::StatisticsType
    Moments1D<TPoint, TStatistics>::
    Dist1(const VectorType* aSumOfPowers, unsigned int aN) const {
        return Dist1(&(this->m_SumOfPowers), this->GetSize(),aSumOfPowers, aN);
    }

    template <class TPoint, class TStatistics>
    typename Moments1D<TPoint, TStatistics>::StatisticsType
    Moments1D<TPoint, TStatistics>::
    Dist1(
    const VectorType* aSumOfPowers1,
    unsigned int aN1,
    const VectorType* aSumOfPowers2,
    unsigned int aN2) const {
        
        StatisticsType vRet = 0;
        StatisticsType vSigma1 = static_cast<StatisticsType>(sqrt((*aSumOfPowers1)[2] / aN1));
        StatisticsType vSigma2 = static_cast<StatisticsType>(sqrt((*aSumOfPowers2)[2] / aN2));
                // Todo: choose the smaller order of both or make order a template param.

        for(int vK = 2; vK < m_Order+1; vK++) {
            double vRooted_moment_this = pow(fabs((*aSumOfPowers1)[vK] / aN1) , 1.0 / vK) / vSigma1;
            double vRooted_moment_arg = pow(fabs((*aSumOfPowers2)[vK] / aN2) , 1.0 / vK) / vSigma2;
            vRet += static_cast<StatisticsType>(fabs(vRooted_moment_this - vRooted_moment_arg));
        }
        return vRet;
    }

    template <class TPoint, class TStatistics>
    typename Moments1D<TPoint, TStatistics>::StatisticsType
    Moments1D<TPoint, TStatistics>::
    Dist2(const Pointer aM) const {
        return Dist2(&(this->m_SumOfPowers), this->GetSize(),
                &(aM->m_SumOfPowers), aM->GetSize());
    }

    template <class TPoint, class TStatistics>
    typename Moments1D<TPoint, TStatistics>::StatisticsType
    Moments1D<TPoint, TStatistics>::
    Dist2(const VectorType* aSumOfPowers, unsigned int aN) const {
        return Dist2(&(this->m_SumOfPowers), this->GetSize(),
                aSumOfPowers, aN);
    }

    template <class TPoint, class TStatistics>
    typename Moments1D<TPoint, TStatistics>::StatisticsType
    Moments1D<TPoint, TStatistics>::
    Dist2(
    const VectorType* aSumOfPowers1,
    unsigned int aN1,
    const VectorType* aSumOfPowers2,
    unsigned int aN2) const {

        StatisticsType vRet = 0;
        StatisticsType vSigma1 = static_cast<StatisticsType>(sqrt((*aSumOfPowers1)[2] / aN1));
        StatisticsType vSigma2 = static_cast<StatisticsType>(sqrt((*aSumOfPowers2)[2] / aN2));
                // Todo: choose the smaller order of both or make order a template param.

        for(int vK = 2; vK < m_Order+1; vK++) {
//            double vRooted_moment_this = pow(fabs((*aSumOfPowers1)[vK] / aN1) / pow(vSigma1,vK), 1.0 / vK);
//            double vRooted_moment_arg = pow(fabs((*aSumOfPowers2)[vK] / aN2) / pow(vSigma2, vK), 1.0 / vK);

            double vRooted_moment_this = pow(fabs((*aSumOfPowers1)[vK] / aN1) , 1.0 / vK) / vSigma1;
            double vRooted_moment_arg = pow(fabs((*aSumOfPowers2)[vK] / aN2) , 1.0 / vK) / vSigma2;
            
            vRet += (vRooted_moment_this - vRooted_moment_arg) *
                    (vRooted_moment_this - vRooted_moment_arg);
        }
        
        vRet = static_cast<StatisticsType>(sqrt(vRet));
        return vRet;
    }

    

    template <class TPoint, class TStatistics>
    typename Moments1D<TPoint, TStatistics>::PointType 
    Moments1D<TPoint, TStatistics>::GetCenterOfMass() {
//        if(m_modified) {
//            m_CenterOfMass = CalcCenterOfMass(m_Points);
//        }
        return m_CenterOfMass;
    }
    
    template <class TPoint, class TStatistics>
    typename Moments1D<TPoint, TStatistics>::PointType 
    Moments1D<TPoint, TStatistics>::GetReferencePoint() {
//        if(m_modified) {
//            m_Centroid = CalcCentroid(m_Points, GetCenterOfMass());
//        }
        return m_CenterOfMass;
    }
    
    template <class TPoint, class TStatistics>
    typename Moments1D<TPoint, TStatistics>::StatisticsType 
    Moments1D<TPoint, TStatistics>::GetMaxDist() {
//        if(m_modified) {
//            Update();
//        }
        return m_MaxSample;
    }
        
    
    template <class TPoint, class TStatistics>
    void Moments1D<TPoint, TStatistics>::
    Update() {

            m_CenterOfMass = CalcCenterOfMass(m_Points);
            
            m_Centroid = CalcCentroid(m_Points, m_CenterOfMass);

            m_Mu = CalcMu(m_Points, m_Centroid);

            /// update the sum of powers and the extremum sample
            m_MaxSample = CalcSumOfPowers(m_Points, m_Centroid, m_Mu, &m_SumOfPowers);

    }

    template <class TPoint, class TStatistics>
    void Moments1D<TPoint, TStatistics>::
    UpdateSumOfPowersWithoutData(
    const VectorType* aSumOfPowers,
    VectorType* aSumOfPowersUpdated,
    unsigned int vNoriginal, double aMuOriginal, double aOldDist) const {

        /// similarly to the update with the point with 2 differences:
        /// 1. we first have to update the mean
        /// 2. The sum of powers is dependent of the lower order sums

        double vNewMean = (aMuOriginal * vNoriginal - aOldDist) / (vNoriginal - 1);
        double vDelta = aOldDist - vNewMean;


        for (unsigned int vP = 0; vP < m_Order + 1; vP++) {
            (*aSumOfPowersUpdated)[vP] = (*aSumOfPowers)[vP];
        }

        /// here we update sequentially since we need the lower order values:
        for (unsigned int vP = 2; vP < m_Order + 1; vP++) { // !! p is the order here
            double vSum = 0;
            for (unsigned int vK = 1; vK <= vP - 2; vK++) { // !! k is the running variable
                vSum += (m_NChooseK[vP][vK]) * (*aSumOfPowersUpdated)[vP - vK] * pow(-vDelta / vNoriginal, vK);
            }
            (*aSumOfPowersUpdated)[vP] -= vSum;

            (*aSumOfPowersUpdated)[vP] -=
                    pow((static_cast<double> (vNoriginal) - 1.0) * vDelta / vNoriginal, vP) *
                    (1.0 - pow(-1.0 / (vNoriginal - 1), static_cast<double> (vP) - 1.0));
        }
    }

    template <class TPoint, class TStatistics>
    void Moments1D<TPoint, TStatistics>::
    UpdateSumOfPowersWithData(
    const VectorType* aSumOfPowers,
    VectorType* aSumOfPowersUpdated,
    unsigned int vN_orig, double aMu_orig, double aNewDistance) const {
        
        unsigned int vN = vN_orig + 1;
        double vDelta = aNewDistance - aMu_orig;
        
        for (unsigned int vP = 0; vP < m_Order + 1; vP++) {
            (*aSumOfPowersUpdated)[vP] = (*aSumOfPowers)[vP];
        }

        for (unsigned int vP = 2; vP < m_Order + 1; vP++) { // !! p is the order here
            double vSum = 0;
            for (unsigned int vK = 1; vK <= vP-2; vK++) { // !! k is the running variable
                vSum += (m_NChooseK[vP][vK]) * (*aSumOfPowers)[vP - vK] * pow(-vDelta / vN, vK);
            }
            (*aSumOfPowersUpdated)[vP] += vSum;
        }

        for (unsigned int vP = 0; vP < m_Order + 1; vP++) { // !! p is the order here
            (*aSumOfPowersUpdated)[vP] += pow((static_cast<double>(vN) - 1.0) * vDelta / vN, vP) *
                    (1.0 - pow(-1.0 / (vN - 1), static_cast<double> (vP) - 1.0));
        }
    }

    template <class TPoint, class TStatistics>
    typename Moments1D<TPoint, TStatistics>::StatisticsType 
    Moments1D<TPoint, TStatistics>::CalcSumOfPowers(
    const PointContainerType* aPoints,
    PointType aReferencePoint, 
    StatisticsType aMu, 
    VectorType* aSumOfPowers) const {
        float vMaxDist = 0;
        for (int vK = 0; vK < m_Order + 1; vK++) {
            (*aSumOfPowers)[vK] = 0;
        }
        /// Get the mean Point (center of mass)
        PointContainerIteratorType vIt = aPoints->begin();
        PointContainerIteratorType vItEnd = aPoints->end();
        for (; vIt != vItEnd; ++vIt) {
            StatisticsType vDist = 0;
            for (unsigned int vD = 0; vD < PointType::Dimension; vD++) {
                vDist += ((*vIt)[vD] - aReferencePoint[vD]) *
                        ((*vIt)[vD] - aReferencePoint[vD]);
            }
            vDist = static_cast<StatisticsType> (sqrt(vDist));

            for (int vK = 2; vK < m_Order + 1; vK++) {
                (*aSumOfPowers)[vK] += static_cast<StatisticsType>(pow(vDist - aMu, vK));
            }
            
            if(vDist > vMaxDist) vMaxDist = vDist;
        }
        return vMaxDist;
    }

    template <class TPoint, class TStatistics>
    typename Moments1D<TPoint, TStatistics>::PointType
    Moments1D<TPoint, TStatistics>::
    CalcCenterOfMass(const PointContainerType* aPoints) const {

        PointType vCOM;
        unsigned int vN = 0;
        /// Get the mean Point (center of mass)
        for (unsigned int vD = 0; vD < PointType::Dimension; vD++) {
            vCOM[vD] = 0;
        }

        /// Get the mean Point (center of mass)
        PointContainerIteratorType vIt = aPoints->begin();
        PointContainerIteratorType vItEnd = aPoints->end();
        for (; vIt != vItEnd; ++vIt) {
            for (unsigned int vD = 0; vD < PointType::Dimension; vD++) {
                vCOM[vD] += (*vIt)[vD];
            }
            vN++;
        }

        for (unsigned int vD = 0; vD < PointType::Dimension; vD++) {
            vCOM[vD] /= vN;
            
        }
        return vCOM;
    }
    
    template <class TPoint, class TStatistics>
    typename Moments1D<TPoint, TStatistics>::PointType
    Moments1D<TPoint, TStatistics>::
    CalcCentroid(const PointContainerType* aPoints, 
            PointType aCenterOfMass) const {
        
        if(m_CentroidType == 0) {
            return aCenterOfMass;
        }
        
        /// Calculate the inertia tensor
         typedef float RealType;
         typedef vnl_matrix<RealType> MatrixType;
         
         MatrixType vTensor(PointType::Dimension,PointType::Dimension,0);
         
         PointContainerIteratorType vIt = (*aPoints).begin();
         PointContainerIteratorType vItEnd = (*aPoints).end();
         if(PointType::Dimension == 2) {
             for (; vIt != vItEnd; ++vIt) {
                 float vX = (*vIt)[0] - aCenterOfMass[0];
                 float vY = (*vIt)[1] - aCenterOfMass[1];
                 
                 vTensor[0][0] += vY*vY;
                 vTensor[0][1] += -vX*vY;
                 vTensor[1][1] += vX*vX;
             }
         } else if (PointType::Dimension == 3) {
             for (; vIt != vItEnd; ++vIt) {
                 float vX = (*vIt)[0] - aCenterOfMass[0];
                 float vY = (*vIt)[1] - aCenterOfMass[1];
                 float vZ = (*vIt)[2] - aCenterOfMass[2];
                 
                 vTensor[0][0] += vY*vY+vZ*vZ;
                 vTensor[0][1] += -vX*vY;
                 vTensor[0][2] += -vX*vZ;
                 vTensor[1][1] += vX*vX+vZ*vZ;
                 vTensor[1][2] += -vY*vZ;
                 vTensor[2][2] += vX*vX+vY*vY;
             }
         }
         /// Symmetrify:
         for (unsigned int vI = 0; vI < PointType::Dimension; vI++) {
             for (unsigned int vJ = 0; vJ < vI; vJ++) {
                 vTensor[vI][vJ] = vTensor[vJ][vI];
             }
         }

         /// normalize the tensor (with uniform weights)
         if(aPoints->size() > 0) vTensor /= aPoints->size();

         /// Solve the eigen-system
         typedef vnl_symmetric_eigensystem<RealType> EigSysType;
         EigSysType vEigSys(vTensor); 
         
         if(m_CentroidType <= PointType::Dimension) {
             /// Get the normalized eigenvector
             vnl_vector<float> vEigV = vEigSys.get_eigenvector(m_CentroidType-1);
             /// scale it according with the sqrt of the corresp. eigenvalue
             vEigV *= sqrt(vEigSys.get_eigenvalue(m_CentroidType-1));
             /// Add the offset vector to the center of mass and return
             PointType vRet; 
             for(unsigned int vD = 0; vD < PointType::Dimension; vD++) {
                 vRet[vD] = aCenterOfMass[vD] + vEigV[vD];
             }
             return vRet;
         } else {
             itk::InvalidArgumentError vE(__FILE__, __LINE__);
             vE.SetLocation("Moments1D");
             vE.SetDescription("Centroid type are valid only from 0--Dimension.");
             throw vE;
         }
         return aCenterOfMass;         
    }

    template <class TPoint, class TStatistics>
    typename Moments1D<TPoint, TStatistics>::StatisticsType
    Moments1D<TPoint, TStatistics>::
    CalcMu(const PointContainerType* aPoints, PointType aReferencePoint) const{
        /// find mu (of the Eucledian distances)
        double vSum = 0;
        double vN = 0;
        PointContainerIteratorType vIt = (*aPoints).begin();
        PointContainerIteratorType vItEnd = (*aPoints).end();
        for (; vIt != vItEnd; ++vIt) {
            StatisticsType vDist = 0;
            for (unsigned int vD = 0; vD < PointType::Dimension; vD++) {
                vDist += ((*vIt)[vD] - aReferencePoint[vD]) *
                        ((*vIt)[vD] - aReferencePoint[vD]);
            }
            vSum += static_cast<StatisticsType> (sqrt(vDist));
            vN++;
        }
        
        return static_cast<StatisticsType>(vSum / vN);
    }


    
    template <class TPoint, class TStatistics>
    void Moments1D<TPoint, TStatistics>::initializeNchooseK(unsigned int aOrder) {
        unsigned int vOrderPlus = aOrder+1;

        Self::m_NChooseK.set_size(vOrderPlus,vOrderPlus);

        //Coefficients
        Self::m_NChooseK.fill(0);
        Self::m_NChooseK(0, 0) = 1;
        Self::m_NChooseK(1, 0) = 1;
        Self::m_NChooseK(1, 1) = 1;
        for (unsigned int i = 2; i < vOrderPlus; i++) {
            Self::m_NChooseK(i, 0) = 1;
            for (unsigned int j = 1; j < i; j++) {
                Self::m_NChooseK(i, j) = Self::m_NChooseK(i - 1, j - 1) + Self::m_NChooseK(i - 1, j);
            }
            Self::m_NChooseK(i, i) = 1;
        }
    }

} // end namespace itk

#endif //_ITKMOMENTS1D_HXX