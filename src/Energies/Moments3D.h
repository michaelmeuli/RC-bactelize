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

#ifndef _MOMENTS3D_H
#define	_MOMENTS3D_H

#include"vnl/vnl_matrix.h"

#include"ShapePrior.h"

#include "LUMatrix.h"

class Moments3D: public ShapePrior
{
private:
        //reference moments
	LUMatrix m_Reference;

        //current mometns and temp moments used by update
	LUMatrix m_Moments;
	LUMatrix m_OnePointMoments;
        LUMatrix m_UpdatedCentralMoments;
        LUMatrix m_UpdatedNormCentralMoments;
        LUMatrix m_UpdatedMoments;
        LUMatrix m_AddPointsMoments;
        LUMatrix m_DeletePointsMoments;
        LUMatrix m_CentralMoments;
        LUMatrix m_NCentralMoments;
        LUMatrix m_TempForCentral;
        vnl_vector<double> m_TempForNorm;
        vnl_vector<double> m_TempForSin;
        vnl_vector<double> m_TempForCos;
        LUMatrix m_dxdydz;
        vnl_vector<double> m_dx;
        vnl_vector<double> m_dy;
        vnl_vector<double> m_dz;
        vnl_matrix<double> m_NChooseK;
        LUMatrix m_GaussianWeight;

        //coefficient of polynomials
        vnl_matrix<double> m_Co_Legendre;
        LUMatrix m_TempForLegendre;


        //for discrete moments
        vnl_matrix<double> m_Co_Chebyshev;
        vnl_matrix<double> m_Table_Chebyshev;
        vnl_matrix<double> m_Table_Legendre;

        //one point change, corresponding measure
	double m_MomentsMeasure;
        //current distance of moments and prior moments
        double m_CurrentMeasure;
public:
	/************************* reference operation *************************/
	void SetReference(LUMatrix *aRef);
	LUMatrix * GetReference();
        /************************* obtain current moments *************************/
        LUMatrix * GetMoments();
        void SetCurrentMeasure(double aM);
        double GetCurrentMeasure();

        /************************* common operations *************************/
        void Init();
        void NChooseK();
        void Co_Legendre();
        void CopyMoments(LUMatrix *src, LUMatrix *dst);
        void NewCenter(double aX, double aY, double aZ, bool AddOrDelete, double *nCenterX, double * nCenterY, double *nCenterZ);
        double AddOnePoint(unsigned int aX, unsigned int aY, unsigned int aZ);
        double DeleteOnePoint(unsigned int aX, unsigned int aY, unsigned int aZ);
        void TranslateMoments(LUMatrix *src, LUMatrix * dst, double dx, double dy, double dz);
        void NormalizeMoments(LUMatrix * src, LUMatrix * dst);
        void ScaleMoments(LUMatrix *src, LUMatrix *dst);
        void ScaleMoments(LUMatrix *src, LUMatrix *dst, double rate);
        void RotateMoments(LUMatrix * src, LUMatrix * dst, double aAngle);
        void LegendreMoments(LUMatrix *GeoMoments, LUMatrix *Legendre);
        void OnePointMoments(double aX, double aY);
        void AddOrDeleteMoments(LUMatrix *Moments, LUMatrix *OnePoint, bool AddOrDelete);
        void GaussianWeight(double aSigma);
        double Distance(LUMatrix *mom1,LUMatrix *mom2, bool WeightedOrNot);
        void Legendre();

        void StepUpdate(vnl_matrix<unsigned int> *aAdd, vnl_matrix<unsigned int> *aDelete);
        void PointsMoments(vnl_matrix<unsigned int> *aPoints, bool AddOrDelete);
	/************************* for continuous moments *************************/
        void CentralMomentsForC();
        void PointPosition(unsigned int aX, unsigned int aY, unsigned int aZ, double *aNX, double *aNY, double *aNZ);
        /************************* for discrete moments *************************/
        void CentralMomentsForD();
        void Co_Chebyshev();
        void Chebyshev();

};

#ifndef ITK_MANUAL_INSTANTIATION
#include "Moments3D.hxx"
#endif
	
#endif	
	
