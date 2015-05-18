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
 * File:   itkMoments1DEnergy.h
 * Author: janickc
 *
 * Created on July 3, 2012, 11:01 AM
 */

#ifndef __ITKMOMENTS1DENERGY_H 
#define	__ITKMOMENTS1DENERGY_H 

#include <set>
#include <vector>
#include "itkObject.h"
#include "itkOrderedPoint.h"
#include "itkMoments1D.h"

namespace itk {
    /// Note: itkMesh resp. itkPointSets do not allow to remove points (without
    ///       knowing a point ID) just by using their coordinates.
    /// Todo: sould probably be inherited from itkDataObject. What functions 
    ///       should be overloaded?
     
    template<class TOrderedPoint, typename TStatistics = double>
    class ITK_EXPORT Moments1DEnergy : public Object{
        
        public:
            
            typedef Moments1DEnergy<TOrderedPoint, TStatistics> Self;
            typedef SmartPointer< Self > Pointer;
            typedef SmartPointer< const Self > ConstPointer;
            
            itkNewMacro(Self);
            
            typedef TOrderedPoint PointType;
            typedef TStatistics StatisticsType;
            typedef itk::Index<PointType::PointDimension> IndexType; 
            

            itkGetConstMacro(NbCentroids, unsigned int);
            itkSetMacro(NbCentroids, unsigned int);
            
            
            void Update();
            
            /** Returns the number of points considered for this moment distribution. */
            unsigned int GetSize() {return m_Points.size();};
            
            private:
                typedef Moments1D<PointType, StatisticsType> Moments1DType;
                typedef typename Moments1DType::Pointer Moments1DPointerType;
        
                typedef std::vector<Moments1DPointerType> Moment1DVectorType;
                Moment1DVectorType m_MomentsObjVec;

                
                typedef std::set<PointType> PointContainerType;
                PointContainerType m_Points;
                
                unsigned int m_NbCentroids;
                unsigned int m_Order;
                std::vector<float> m_Coeff;                
                bool m_modified;
                
                public:
                    void Init(unsigned int aOrder, unsigned int aNbCentroids);
                    void AddPoint(PointType);
                    void RemovePoint(PointType);
                    void AddPoint(IndexType);
                    void RemovePoint(IndexType);
                    void SetCoefficient(unsigned int aIndex, float aCoeff);    
                    
                    StatisticsType Dist1(Pointer aMomentsEnergy);
                    StatisticsType EnergyDiffL1Approx(const Pointer aMomentsEnergy, PointType aPoint, bool aAdd);
                    

    }; // end class
} // end namespace

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMoments1DEnergy.hxx"
#endif

#endif	/* __ITKMOMENTS1DENERGY_H */

