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
/* 
 * File:   itkMoments1D.h
 * Author: janickc
 *
 * Created on October 12, 2011, 3:11 PM
 */

#ifndef _ITKMOMENTS1D_H
#define	_ITKMOMENTS1D_H

#include "itkImage.h"
#include "vnl/vnl_matrix.h"
#include <set>

namespace itk {

    // It is unclear how to identify points when using ITKs PointSet using only
    // a 1D unsigned integer (called IdentifierType). So we use
    // a std::set here (where keys are values at the same time) instead to 
    // store the points.
    // There is no ITK wrapper for the std::set (there is one for std::map).

    // TODO should most probably a itkDataObject OR part of the statistic module OR
    //      a PointSetMetric OR a SinglevaluedCostFunction.

    // TODO: point type should actually be of type integer such that they can be
    //       later retrieved in the set. This means that we need an internal
    //       point type of with TCoord of statisticstype.

    // TODO: the update shouldn't be invoked by the user. It should be of private scope.
    //       The update should rather be invoked only when a distance is requested AND
    //       the points indeed had been changed. This needs a modified flag.

    // TODO: the dist measurements should be ind a different class (derived from
    //       itkMetric) wheras this class should be some specialized point container.
    //       The output of this class should be the vector of normalized moments.

    // TODO: To not store the points, maybe make an method called 'aggregate' which
    //       takes a point as input and accumulates all the means (COM, mu, N) etc.
    //       One could potentially use boosts aggregates.



    template<class TPoint, typename TStatistics = double>
    class ITK_EXPORT Moments1D : public Object{
    public:
//        typedef Moments1D Self;
        typedef Moments1D<TPoint, TStatistics> Self;
//        typedef Object Superclass;
        typedef SmartPointer< Self > Pointer;
        typedef SmartPointer< const Self > ConstPointer;

        itkNewMacro(Self);

        Moments1D(){};

        typedef TStatistics StatisticsType;
        typedef itk::Index<TPoint::PointDimension> IndexType; 
//        typedef TPointSet PointSetType;
//        typedef typename PointSetType::PointType PointType;
//        typedef typename PointSetType::Pointer PointSetPointerType;

        typedef TPoint PointType;
        typedef std::set<PointType> PointContainerType;
        typedef PointContainerType* PointContainerPointerType;
        typedef typename PointContainerType::iterator PointContainerIteratorType;
        
        PointContainerPointerType m_Points;
        itkSetMacro(Points, PointContainerPointerType);

        typedef vnl_vector<StatisticsType> VectorType;
        typedef vnl_matrix<StatisticsType> MatrixType;
        
//        void operator+=(typename PointSetType::PointType);
//        void operator-=(typename PointSetType::PointType);

        
        void Update();
        void Init(unsigned int aOrder);
        void Init(unsigned int aOrder, unsigned char aType);
        

//        itkGetConstMacro(N, unsigned int);
//        itkSetMacro(Order, unsigned int);
//        itkGetConstMacro(Order, unsigned int);
        itkGetConstReferenceMacro(SumOfPowers, VectorType);

        unsigned int GetSize() const;

        StatisticsType Dist1(const Pointer aM) const;
        StatisticsType Dist1(const VectorType* aSumOfPowers, unsigned int aN) const;
        StatisticsType Dist1(
                const VectorType* aSumOfPowers1,
                unsigned int aN1,
                const VectorType* aSumOfPowers2,
                unsigned int aN2) const;
        StatisticsType Dist2(const Pointer aM) const;
        StatisticsType Dist2(const VectorType* aSumOfPowers, unsigned int aN) const;
        StatisticsType Dist2(
                const VectorType* aSumOfPowers1,
                unsigned int aN1,
                const VectorType* aSumOfPowers2,
                unsigned int aN2) const;


        void CalcSumOfPowersWithPoint(PointType aPoint, VectorType* aSumOfPowersOutput);
        void CalcSumOfPowersWithoutPoint(PointType aPoint, VectorType* aSumOfPowersOutput);
        void CalcSumOfPowersApproxWithPoint(PointType aPoint, VectorType* aSumOfPowersOutput) const;
        void CalcSumOfPowersApproxWithoutPoint(PointType aPoint, VectorType* aSumOfPowersOutput) const;

        
        PointType GetCenterOfMass();
        PointType GetReferencePoint();
        StatisticsType GetMaxDist();
    private:
        StatisticsType m_Mu;
        VectorType m_SumOfPowers; // from here the moments can be calculated
        StatisticsType m_MaxSample;
        
        unsigned int m_Order;

        PointType m_CenterOfMass;
        PointType m_Centroid;
        unsigned char m_CentroidType;
        static bool m_NchooseK_initialized;
        static MatrixType m_NChooseK;

        void initializeNchooseK(unsigned int aOrder);
        PointType CalcCenterOfMass(const PointContainerType* aPoints) const;
        StatisticsType CalcMu(const PointContainerType* aPoints, PointType aReferencePoint) const;
        PointType CalcCentroid(const PointContainerType* aPoints, PointType aCenterOfMass) const;

        void UpdateSumOfPowersWithData(
                const VectorType* aSumOfPowers,
                VectorType* aSumOfPowersUpdated,
                unsigned int vNoriginal,
                double aMuOriginal,
                double aNewDist) const;
        void UpdateSumOfPowersWithoutData(
                const VectorType* aSumOfPowers,
                VectorType* aSumOfPowersUpdated,
                unsigned int vNoriginal,
                double aMuOriginal,
                double aOldDist) const;
        StatisticsType CalcSumOfPowers(
                const PointContainerType * aPoints,
                PointType aReferencePoint,
                StatisticsType aMu,
                VectorType* aSumOfPowers) const;

    };

template<class TPoint, typename TStatistics>
ITK_EXPORT std::ostream & operator<<(std::ostream & os,
                                     const Moments1D< TPoint, TStatistics > & aMoments) {
        for (unsigned int vK = 0; vK < aMoments.getOrder() + 1; aMoments) {
            os << "SumOfPowers vecotr:\n " << aMoments;
            os << std::endl;
        }
    return os;
}

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMoments1D.hxx"
#endif

#endif	/* ITKMOMENTS1D_H */

