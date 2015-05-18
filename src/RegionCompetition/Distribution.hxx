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
 *  Distribution.hxx
 *  HelloWorld
 *
 *  Created by Yuanhao Gong on 9/17/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __Distribution_hxx
#define __Distribution_hxx
#include<iostream>
#include "vnl/vnl_cross.h"

void Distribution::SetNumberOfBins(int aS)
{
	m_NumberOfBins=aS;
}

int Distribution::GetNumberOfBins()
{
	return m_NumberOfBins;
}

void Distribution::SetRepeat(int aS)
{
	m_Repeat=aS;
}

int Distribution::GetRepeat()
{
	return m_Repeat;
}

void Distribution::SetSampeSize(int aS)
{
	m_SampleSize=aS;
}

int Distribution::GetSampleSize()
{
	return m_SampleSize;
}

void Distribution::SetReference(vnl_vector<double> * aRef)
{
	m_Reference=*aRef;
}

vnl_vector<double> * Distribution::GetReference()
{
	return &m_Reference;
}

vnl_vector<double>* Distribution::GetA3()
{
	return &m_A3;
}

vnl_vector<double>* Distribution::GetD1()
{
	return &m_D1;
}

vnl_vector<double>* Distribution::GetD2()
{
	return &m_D2;
}

vnl_vector<double>* Distribution::GetD3()
{
	return &m_D3;
}

vnl_vector<double>* Distribution::GetD4()
{
	return &m_D4;
}

/************************* operation functions *************************/
void Distribution::Initial()
{
	m_NumberOfBins=20;
	m_Repeat=4;
	m_SampleSize=GetSize();
	
	m_PDF.set_size(m_NumberOfBins);
	m_Reference.set_size(m_NumberOfBins);
	m_A3.set_size(m_NumberOfBins);
	m_D1.set_size(m_NumberOfBins);
	m_D2.set_size(m_NumberOfBins);
	m_D3.set_size(m_NumberOfBins);
	m_D4.set_size(m_NumberOfBins);
	m_Container.set_size(m_NumberOfBins, m_Repeat);
	m_Temp.set_size(m_NumberOfBins);
}

void Distribution::Initial(int aNb, int aRepeat, int aSampleSize)
{
	m_NumberOfBins=aNb;
	m_Repeat=aRepeat;
	m_SampleSize=aSampleSize;
	
	m_PDF.set_size(m_NumberOfBins);
	m_Reference.set_size(m_NumberOfBins);
	m_A3.set_size(m_NumberOfBins);
	m_D1.set_size(m_NumberOfBins);
	m_D2.set_size(m_NumberOfBins);
	m_D3.set_size(m_NumberOfBins);
	m_D4.set_size(m_NumberOfBins);
	m_Container.set_size(m_NumberOfBins, m_Repeat);
}

void Distribution::ComputeA3()
{
	const int vS=GetSize();
	const int vD=GetDimension();
	double vStep=ShapePrior::PI/(m_NumberOfBins-1);
	
	//inner temp var
	double temp=0;
	unsigned int a,b,c;
	double d1, d2;
	vnl_vector<double> V1(vD), V2(vD);
	int index;
	
	//container
	m_Container.fill(0);
	
	//random number generator
	vnl_random w(vS);
	
	for (int i=0; i<m_Repeat; i++) 
	{
		//vcl_cerr<<i<<"**********************************************"<<std::endl;

		for (int j=0; j<m_SampleSize; j++)
		{
			a=w.lrand32(vS-1);
			b=w.lrand32(vS-1);
			c=w.lrand32(vS-1);
			
			//vcl_cerr<<j<<" "<<a<<" "<<b<<" "<<c<<std::endl;
			V1=(*GetInputData()).get_row(b)-(*GetInputData()).get_row(a);
			V2=(*GetInputData()).get_row(b)-(*GetInputData()).get_row(c);
			d1=V1.two_norm();
			d2=V2.two_norm();
			//remove numerical unstable points
			if ((d1>ShapePrior::MAX_ZERO_ERROR) && (d2>ShapePrior::MAX_ZERO_ERROR) ) 
			{
				temp=dot_product(V1, V2)/(d1*d2);
				//vcl_cerr<<j<<" dot product "<<temp<<std::endl;
				
				//deal with interval boundery, due to numerical reasons
				if (fabs(temp-1)<ShapePrior::MAX_ZERO_ERROR) {
					temp=0;
				}else {
					if (fabs(temp+1)<ShapePrior::MAX_ZERO_ERROR) {
						temp=ShapePrior::PI/2;
					}else {
						temp=acos(temp);
					}

				}				
				//vcl_cerr<<j<<" "<<temp<<std::endl;
				index=round(temp/vStep);
				if (index>=m_NumberOfBins || index <0) {
					std::cout<<"Index error in A3, please check."<<std::endl;
				}
				//vcl_cerr<<j<<" "<<index<<std::endl;
				m_Container(index,i)=m_Container(index,i)+1;
			}
		}
	}
	
	//get everage of repeat
	m_A3.fill(0);
	for (int i=0; i<m_Repeat; i++)
	{
		m_A3=m_A3+m_Container.get_column(i)/m_Container.get_column(i).one_norm();
	}
	
	double dTemp=m_Repeat*1.0;
	m_A3/=dTemp;
}

void Distribution::ComputeD1()
{
	const int vS=GetSize();
	const int vD=GetDimension();
	
	//translate to the center
	vnl_matrix<double> vM(vS,vD);
	for (int i=0; i<vD; i++) {
		vM.set_column(i,(*GetInputData()).get_column(i)-(*GetInputData()).get_column(i).mean());
	}
	//distance to center
	vnl_vector<double> vT(vS);
	
	vT.fill(0);
	
	for (int i=0; i<vS; i++) {
		vT=vT+element_product(vM.get_column(i), vM.get_column(i));
	}
	//max distance
	double vMax=vT.max_value();
	vT/=vMax;
	
	//step
	double step=1.0/(m_NumberOfBins-1);
	
	//container
	vnl_vector<double> C(m_NumberOfBins);
	C.fill(0);
	int index=0;
	for (int i=0; i<vS; i++) {
		index=round(vT(i)/step);
		C(index)=C(index)+1;
	}
	m_D1=C/C.one_norm();
}

void Distribution::ComputeD2()
{
	const int vS=GetSize();
	const int vD=GetDimension();
	double vStep=1.0/(m_NumberOfBins-1.0);
	
	//inner temp var
	unsigned int a,b;
	double d1;
	vnl_vector<double> V1(vD);
	int index;
	
	//container
	m_Container.fill(0);
	
	//random number generator
	vnl_random w(vS);
	
	vnl_vector<double> vContainVector(m_SampleSize);
	for (int i=0; i<m_Repeat; i++) 
	{
		//vcl_cerr<<i<<"**********************************************"<<std::endl;

		for (int j=0; j<m_SampleSize; j++)
		{
			a=w.lrand32(vS-1);
			b=w.lrand32(vS-1);
			
			//vcl_cerr<<j<<" "<<a<<" "<<b<<" "<<c<<std::endl;
			V1=(*GetInputData()).get_row(b)-(*GetInputData()).get_row(a);
			d1=V1.two_norm();
			vContainVector(j)=d1;
		}
		vContainVector/=vContainVector.max_value();
		for(int j=0;j<m_SampleSize;j++)
		{
			//vcl_cerr<<j<<" "<<temp<<std::endl;
			index=round(vContainVector(j)/vStep);
			if (index>=m_NumberOfBins || index <0) {
				std::cout<<"Index error in D2, please check."<<std::endl;
			}
			//vcl_cerr<<j<<" "<<index<<std::endl;
			m_Container(index,i)=m_Container(index,i)+1;
		}
	}
	
	//get everage of repeat
	m_D2.fill(0);
	for (int i=0; i<m_Repeat; i++)
	{
		m_D2=m_D2+m_Container.get_column(i)/m_Container.get_column(i).one_norm();
	}
	
	m_D2/=(m_Repeat*1.0);	
}

void Distribution::ComputeD3()
{
	const int vS=GetSize();
	const int vD=GetDimension();
	double vStep=1.0/(m_NumberOfBins-1);
	
	//inner temp var
	double temp=0;
	unsigned int a,b,c;
	double d1;
	vnl_vector<double> V1(vD), V2(vD);
	int index;
	
	//container
	m_Container.fill(0);
	
	//random number generator
	vnl_random w(vS);
	vnl_vector<double> vContainVector(m_SampleSize);
	
	for (int i=0; i<m_Repeat; i++) 
	{
		//vcl_cerr<<i<<"**********************************************"<<std::endl;

		for (int j=0; j<m_SampleSize; j++)
		{
			a=w.lrand32(vS-1);
			b=w.lrand32(vS-1);
			c=w.lrand32(vS-1);
			
			//vcl_cerr<<j<<" "<<a<<" "<<b<<" "<<c<<std::endl;
			V1=(*GetInputData()).get_row(b)-(*GetInputData()).get_row(a);
			V2=(*GetInputData()).get_row(b)-(*GetInputData()).get_row(c);
                        if ((*GetInputData()).columns()==2)
                        {
                            vContainVector(j)=fabs(vnl_cross_2d(V1,V2));
                        }else{
                            vContainVector(j)=fabs(vnl_cross_3d(V1,V2).magnitude());
                        }
		}
		vContainVector/=vContainVector.max_value();
		for(int j=0;j<m_SampleSize;j++)
		{				
			temp=vContainVector(j);
			//vcl_cerr<<j<<" "<<temp<<std::endl;
			index=round(temp/vStep);
			if (index>=m_NumberOfBins || index <0) 
			{
				std::cout<<"Index error in D3, please check."<<std::endl;
			}
			//vcl_cerr<<j<<" "<<index<<std::endl;
			m_Container(index,i)=m_Container(index,i)+1;
		}
	}
	
	//get everage of repeat
	m_D3.fill(0);
	for (int i=0; i<m_Repeat; i++)
	{
		m_D3=m_D3+m_Container.get_column(i)/m_Container.get_column(i).one_norm();
	}
	
	d1=m_Repeat*1.0;
	m_D3/=d1;

}

void Distribution::ComputeD4()
{
	const int vS=GetSize();
	const int vD=GetDimension();
	double vStep=1.0/(m_NumberOfBins-1);
	
	//inner temp var
	double temp=0;
	unsigned int a,b,c,d;
	double d1;
	vnl_vector<double> V1(vD), V2(vD), V3(vD);
	int index;
	
	//container
	m_Container.fill(0);
	
	vnl_matrix_fixed<double,3,3> vInnerMatrix;
	//random number generator
	vnl_random w(vS);
	vnl_vector<double> vContainVector(m_SampleSize);	
	for (int i=0; i<m_Repeat; i++) 
	{
		//vcl_cerr<<i<<"**********************************************"<<std::endl;
		for (int j=0; j<m_SampleSize; j++)
		{
			a=w.lrand32(vS-1);
			b=w.lrand32(vS-1);
			c=w.lrand32(vS-1);
			d=w.lrand32(vS-1);
			
			//vcl_cerr<<j<<" "<<a<<" "<<b<<" "<<c<<std::endl;
			V1=(*GetInputData()).get_row(b)-(*GetInputData()).get_row(a);
			V2=(*GetInputData()).get_row(c)-(*GetInputData()).get_row(a);
			V3=(*GetInputData()).get_row(d)-(*GetInputData()).get_row(a);
			vInnerMatrix.set_row(0,V1);
			vInnerMatrix.set_row(1,V2);
			vInnerMatrix.set_row(2,V3);
			vContainVector(j)=fabs(vnl_det(vInnerMatrix));
		}
		vContainVector/=vContainVector.max_value();
		for(int j=0;j<m_SampleSize;j++)
		{				
			temp=vContainVector(j);
			//vcl_cerr<<j<<" "<<temp<<std::endl;
			index=round(temp/vStep);
			if (index>=m_NumberOfBins || index <0) 
			{
				std::cout<<"Index error in D3, please check."<<std::endl;
			}
			//vcl_cerr<<j<<" "<<index<<std::endl;
			m_Container(index,i)=m_Container(index,i)+1;
		}
	}
	
	//get everage of repeat
	m_D4.fill(0);
	for (int i=0; i<m_Repeat; i++)
	{
		m_D4=m_D4+m_Container.get_column(i)/m_Container.get_column(i).one_norm();
	}
	
	d1=m_Repeat*1.0;
	m_D4/=d1;

}

double Distribution::Distance(vnl_vector<double> * aPDF1, vnl_vector<double> *aPDF2)
{
//	int vS=(*aPDF1).size();
	double temp=0;
	if ((*aPDF1).size()!=(*aPDF2).size())
	{
		std::cout<<"Error in using Distance of Shape Distribution:" <<std::endl;
		return 0;
	}
	else
	{
		m_Temp=(*aPDF1)-(*aPDF2);
		temp=m_Temp.two_norm();
	}
	
	return temp;
}

void Distribution::Compute()
{
	ComputeA3();
	ComputeD1();
	ComputeD2();
	ComputeD3();
}

#endif