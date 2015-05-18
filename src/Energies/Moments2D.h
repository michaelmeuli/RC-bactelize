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
 * Yuanhao Gong
 *=========================================================================*/
#ifndef _MOMENTS2D_H
#define	_MOMENTS2D_H

#include"vnl/vnl_matrix.h"
#include"vnl/vnl_gamma.h"
#include"vnl/algo/vnl_svd_economy.h"
#include"vnl/algo/vnl_svd.h"
#include"ShapePrior.h"

namespace itk {

class Moments2D : public ShapePrior
{
public:
    typedef Moments2D Self;
    typedef itk::SmartPointer< Self > Pointer;
    typedef itk::SmartPointer< const Self > ConstPointer;
    
    itkNewMacro(Self);
    
private:
        //reference moments
	vnl_matrix<double> m_Reference;
        vnl_matrix<double> m_ReferenceNCentralMoments;
        
        //current mometns and temp moments used by update
	vnl_matrix<double> m_Moments;
	vnl_matrix<double> m_OnePointMoments;
        vnl_matrix<double> m_UpdatedCentralMoments;
        vnl_matrix<double> m_UpdatedNormCentralMoments;
        vnl_matrix<double> m_UpdatedMoments;
        vnl_matrix<double> m_AddPointsMoments;
        vnl_matrix<double> m_DeletePointsMoments;
        vnl_matrix<double> m_CentralMoments;
        vnl_matrix<double> m_NCentralMoments;
        vnl_matrix<double> m_RCentralMoments;
        vnl_matrix<double> m_TempForCentral;
        vnl_matrix<double> m_TempForRotation; 
        vnl_vector<double> m_TempForNorm;
        vnl_vector<double> m_TempForSin;
        vnl_vector<double> m_TempForCos;
        vnl_matrix<double> m_dxdy;
        vnl_vector<double> m_dx;
        vnl_vector<double> m_dy;
        vnl_matrix<double> m_NChooseK;
        vnl_matrix<double> m_GaussianWeight;
        
        //coefficient of polynomials
        vnl_matrix<double> m_Co_Legendre;
        vnl_matrix<double> m_TempForLegendre;
        

        //for discrete moments
        vnl_matrix<double> m_Co_Chebyshev;
        vnl_matrix<double> m_Table_Chebyshev;
        vnl_matrix<double> m_Table_Legendre;

        //for angle integration
        vnl_matrix<double> m_Integrate_Angle;
        //one point change, corresponding measure
	double m_MomentsMeasure;
        //current distance of moments and prior moments
        double m_CurrentMeasure;
        //scale ratio for scaled moments
        double m_MomentsScaleRatio;
        //temp for axes rotation matrix
        vnl_matrix_fixed<double,3,3> m_TempRotationMatrix;
public:
	/************************* reference operation *************************/
	void SetReference(vnl_matrix<double> *aRef);
	vnl_matrix<double> * GetReference();
        void ComputeReference();
        /************************* obtain current moments *************************/
        vnl_matrix<double> * GetMoments();
        void SetCurrentMeasure(double aM);
        double GetCurrentMeasure();

        /************************* common operations *************************/
        void Init();
        void NChooseK();
        void Co_Legendre();
        void Co_Integrate_Angle();
        void CopyMoments(vnl_matrix<double> *src, vnl_matrix<double> *dst);
        void NewCenter(double aX, double aY, bool AddOrDelete, double *nCenterX, double * nCenterY);
        double AddOnePoint(unsigned int aX, unsigned int aY);
        double DeleteOnePoint(unsigned int aX, unsigned int aY);
        void TranslateMoments(vnl_matrix<double> *src, vnl_matrix<double> * dst, double dx, double dy);
        void NormalizeMoments(vnl_matrix<double> * src, vnl_matrix<double> * dst);
        void ScaleMoments(vnl_matrix<double> *src, vnl_matrix<double> *dst);
        void ScaleMoments(vnl_matrix<double> *src, vnl_matrix<double> *dst, double rate);
        void RotateMoments(vnl_matrix<double> * src, vnl_matrix<double> * dst, double aAngle);
        void IntegrateAngle(vnl_matrix<double> *src, vnl_matrix<double> *dst);
        void LegendreMoments(vnl_matrix<double> *GeoMoments, vnl_matrix<double> *Legendre);
        void OnePointMoments(double aX, double aY);
        void AddOrDeleteMoments(vnl_matrix<double> *Moments, vnl_matrix<double> *OnePoint, bool AddOrDelete);
        void GaussianWeight(double aSigma);
        double Distance(vnl_matrix<double> *mom1,vnl_matrix<double> *mom2, bool WeightedOrNot);
        void Legendre();
        

        void StepUpdate(vnl_matrix<unsigned int> *aAdd, vnl_matrix<unsigned int> *aDelete);
        void PointsMoments(vnl_matrix<unsigned int> *aPoints, bool AddOrDelete);
	/************************* for continuous moments *************************/
        void CentralMomentsForC();
        void PointPosition(unsigned int aX, unsigned int aY, double *aNX, double *aNY);
        /************************* for discrete moments *************************/
        void CentralMomentsForD();
        void Table_Chebyshev(unsigned int aSize, unsigned int aOrder);
        void ChebyshevMoments(vnl_matrix<unsigned int> * aCoord);
        void Chebyshev();
        
        //estimate the best angle that the moment1 should rotate to match moment2
        double TheBestAngle(vnl_matrix<double> *aMom1, vnl_matrix_fixed<double,3,3> *aVectors1, vnl_matrix<double> *aMom2 , vnl_matrix_fixed<double,3,3> *aVectors2);
        
        //compute the rotation matrix
        vnl_matrix_fixed<double,3,3> RotationMatrixFromMoments(vnl_matrix<double> * CentralMoments);
        
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "Moments2D.hxx"
#endif
	
#endif	
	
