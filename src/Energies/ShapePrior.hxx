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
#ifndef __ShapePrior_hxx
#define __ShapePrior_hxx

namespace itk {

double ShapePrior::PI=3.14159265358979323846;
double ShapePrior::MAX_ZERO_ERROR=0.0000001;

void ShapePrior::SetInputData(vnl_matrix<unsigned int> *aData, unsigned int aWidth, unsigned int aHight, unsigned int aDepth)
{
	m_Coordinate=aData;
	m_DataSize=(*m_Coordinate).rows();
	m_Dimension=(*m_Coordinate).columns();	
	m_HalfWidth=aWidth/2;
	m_HalfHight=aHight/2;
        //the lareger one
        m_ImageSize=aWidth>aHight?aWidth:aHight;
        if (m_Dimension==3)
        {
            m_ImageSize=m_ImageSize>aDepth?m_ImageSize:aDepth;
        }
        //twice imagesize, just make sure that when translation, rotation, scale,
        //the coordinate is till ligal
        m_Scale=2*m_ImageSize;
        m_ScaleDouble=(double)m_Scale;
        unsigned int vx, vy, vz;
        vx=0;
        vy=0;
        vz=0;
        unsigned int ** pData=m_Coordinate->data_array();
        if (m_Dimension==2)
        {
            for (unsigned int i=0;i<m_DataSize;i++)
            {
                vx+=pData[i][0];
                vy+=pData[i][1];
            }
            m_CenterX=vx/m_DataSize;
            m_CenterXDouble=((double)vx+0.0)/((double)m_DataSize);
            m_CenterY=vy/m_DataSize;
            m_CenterYDouble=((double)vy+0.0)/((double)m_DataSize);

        }
        if (m_Dimension==3)
        {
            for (unsigned int i=0;i<m_DataSize;i++)
            {
                vx+=pData[i][0];
                vy+=pData[i][1];
                vz+=pData[i][2];
            }
            m_CenterX=vx/m_DataSize;
            m_CenterXDouble=((double)vx+0.0)/((double)m_DataSize);
            m_CenterY=vy/m_DataSize;
            m_CenterYDouble=((double)vy+0.0)/((double)m_DataSize);
            m_CenterZ=vz/m_DataSize;
            m_CenterZDouble=((double)vz+0.0)/((double)m_DataSize);
        }

}


vnl_matrix<unsigned int> * ShapePrior::GetInputData()
{
	return m_Coordinate;
}

void ShapePrior::SetOrder(unsigned int aOrder)
{
	m_Order=aOrder;
}

unsigned int ShapePrior::GetOrder()
{
	return m_Order;
}

void ShapePrior::SetDimension(unsigned int aDimension)
{
	m_Dimension=aDimension;
}

unsigned int ShapePrior::GetDimension()
{
	return m_Dimension;
}

void ShapePrior::SetDataSize(unsigned int aSize)
{
	m_DataSize=aSize;
}

unsigned int ShapePrior::GetDataSize()
{
	return m_DataSize;
}

void ShapePrior::SetImageSize(unsigned int aSize)
{
	m_ImageSize=aSize;
}

unsigned int ShapePrior::GetImageSize()
{
	return m_ImageSize;
}

void ShapePrior::SetCenterX(unsigned int aX)
{
	m_CenterX=aX;
}

unsigned int ShapePrior::GetCenterX()
{
	return m_CenterX;
}

void ShapePrior::SetCenterY(unsigned int aY)
{
	m_CenterY=aY;
}

unsigned int ShapePrior::GetCenterY()
{
	return m_CenterY;
}

void ShapePrior::SetCenterZ(unsigned int aZ)
{
	m_CenterZ=aZ;
}

unsigned int ShapePrior::GetCenterZ()
{
	return m_CenterZ;
}

void ShapePrior::SetCenterXDouble(double aX)
{
    m_CenterXDouble=aX;
}

double ShapePrior::GetCenterXDouble()
{
    return m_CenterXDouble;
}

void ShapePrior::SetCenterYDouble(double aY)
{
    m_CenterYDouble=aY;
}

double ShapePrior::GetCenterYDouble()
{
    return m_CenterYDouble;
}

void ShapePrior::SetCenterZDouble(double aZ)
{
    m_CenterZDouble=aZ;
}

double ShapePrior::GetCenterZDouble()
{
    return m_CenterZDouble;
}

void ShapePrior::SetScale(unsigned int aS)
{
	m_Scale=aS;
}

unsigned int ShapePrior::GetScale()
{
	return m_Scale;
}

void ShapePrior::SetScaleDouble(double aS)
{
	m_ScaleDouble=aS;
}

double ShapePrior::GetScaleDouble()
{
	return m_ScaleDouble;
}

void ShapePrior::SetTranslationX(unsigned int aX)
{
	m_TranslationX=aX;
}

unsigned int ShapePrior::GetTranslationX()
{
	return m_TranslationX;
}

void ShapePrior::SetTranslationY(unsigned int aY)
{
	m_TranslationY=aY;
}

unsigned int ShapePrior::GetTranslationY()
{
	return m_TranslationY;
}

void ShapePrior::SetTranslationZ(unsigned int aZ)
{
	m_TranslationZ=aZ;
}

unsigned int ShapePrior::GetTranslationZ()
{
	return m_TranslationZ;
}

void ShapePrior::SetTranslationXDouble(double aX)
{
    m_TranslationXDouble=aX;
}

double ShapePrior::GetTranslationXDouble()
{
    return m_TranslationXDouble;
}

void ShapePrior::SetTranslationYDouble(double aY)
{
    m_TranslationYDouble=aY;
}

double ShapePrior::GetTranslationYDouble()
{
    return m_TranslationYDouble;
}

void ShapePrior::SetTranslationZDouble(double aZ)
{
    m_TranslationZDouble=aZ;
}

double ShapePrior::GetTranslationZDouble()
{
    return m_TranslationZDouble;
}

void ShapePrior::SetAngleA(double aA)
{
	m_AngleA=aA;
}

double ShapePrior::GetAngleA()
{
	return m_AngleA;
}

void ShapePrior::SetAngleB(double aB)
{
	m_AngleB=aB;
}

double ShapePrior::GetAngleB()
{
	return m_AngleB;
}
void ShapePrior::SetPriorAxesVector(vnl_matrix_fixed<double,3,3> * aPrior_Axes_Vectors)
{
    for (unsigned int i = 0;i<3;i++)
        for (unsigned int j = 0;j<3;j++)
            m_Prior_Axes_Vectors[i][j] = (*aPrior_Axes_Vectors)[i][j];
}

vnl_matrix_fixed<double,3,3> * ShapePrior::GetPriorAxesVector()
{
    return &m_Prior_Axes_Vectors;
}


} // end namespace itk
#endif
