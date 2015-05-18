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
#ifndef _SHAPEPRIOR_H
#define	_SHAPEPRIOR_H

//std head
#include <cmath>
//vnl head
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_matrix_fixed.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_random.h"
#include "vnl/vnl_det.h"

namespace itk {

class ShapePrior : public itk::Object{
    
public:
    typedef ShapePrior Self;
    typedef SmartPointer< Self > Pointer;
    typedef SmartPointer< const Self > ConstPointer;

    itkNewMacro(Self);
    
private:
	/************************* just record data from image *************************/
	unsigned int m_Order;
	unsigned int m_Dimension;
        unsigned int m_DataSize;
	unsigned int m_ImageSize;//require the image be in squire format

	//original data
	vnl_matrix<unsigned int> * m_Coordinate;
	
	//center of the processed data
	unsigned int m_CenterX;
	unsigned int m_CenterY;
	unsigned int m_CenterZ;
        double m_CenterXDouble;
	double m_CenterYDouble;
        double m_CenterZDouble;

        //scale parameters, int and float scale parameters
	unsigned int m_Scale;
        double m_ScaleDouble;
	
	//translation parameters, just record information from data
	unsigned int m_TranslationX;
	unsigned int m_TranslationY;
	unsigned int m_TranslationZ;
        double m_TranslationXDouble;
        double m_TranslationYDouble;
        double m_TranslationZDouble;

        //rotation angles
        double m_AngleA;
        double m_AngleB;
        //template axes vectors
        vnl_matrix_fixed<double,3,3> m_Prior_Axes_Vectors;
public:	
	/************************* global var *************************/
	static double PI;
	static double MAX_ZERO_ERROR;
        //half width and hight of the image
	unsigned int m_HalfWidth;
	unsigned int m_HalfHight;
	/************************* set and get functions *************************/
	void SetInputData(vnl_matrix<unsigned int> *aData, unsigned int aWidth, unsigned int aHight, unsigned int aDepth);
	vnl_matrix<unsigned int> * GetInputData();
	
	void SetOrder(unsigned int aOrder);
	unsigned int GetOrder();

        void SetDimension(unsigned int aDimension);
        unsigned int GetDimension();

        void SetDataSize(unsigned int aSize);
        unsigned int GetDataSize();

        void SetImageSize(unsigned int aSize);
        unsigned int GetImageSize();
	
	void SetCenterX(unsigned int aX);
	unsigned int GetCenterX();

	void SetCenterY(unsigned int aY);
	unsigned int GetCenterY();

	void SetCenterZ(unsigned int aZ);
	unsigned int GetCenterZ();

        void SetCenterXDouble(double aX);
	double GetCenterXDouble();

        void SetCenterYDouble(double aY);
	double GetCenterYDouble();

        void SetCenterZDouble(double aZ);
	double GetCenterZDouble();

	void SetScale(unsigned int aM);
	unsigned int GetScale();

        void SetScaleDouble(double aM);
        double GetScaleDouble();

	void SetTranslationX(unsigned int aX);
	unsigned int GetTranslationX();

	void SetTranslationY(unsigned int aY);
	unsigned int GetTranslationY();

	void SetTranslationZ(unsigned int aZ);
	unsigned int GetTranslationZ();

        void SetTranslationXDouble(double aX);
	double GetTranslationXDouble();

        void SetTranslationYDouble(double aY);
	double GetTranslationYDouble();

        void SetTranslationZDouble(double aZ);
	double GetTranslationZDouble();

        void SetAngleA(double aAngle);
        double GetAngleA();

        void SetAngleB(double aAngle);
        double GetAngleB();

        void SetPriorAxesVector(vnl_matrix_fixed<double,3,3> * aPrior_Axes_Vectors);
        vnl_matrix_fixed<double,3,3> * GetPriorAxesVector();
	
};
} // end namespace itk
#ifndef ITK_MANUAL_INSTANTIATION
#include "ShapePrior.hxx"
#endif

#endif
