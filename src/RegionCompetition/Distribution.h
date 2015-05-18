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
 *  DISTRIBUTION.h
 *  HelloWorld
 *
 *  Created by Yuanhao Gong on 9/17/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _DISTRIBUTION_H
#define	_DISTRIBUTION_H

#include"vnl/vnl_matrix.h"

#include"ShapePrior.h"

class Distribution: public ShapePrior
{
private:
	int m_NumberOfBins;
	int m_Repeat;
	int m_SampleSize;
	vnl_vector<double> m_PDF;
	vnl_vector<double> m_Reference;
	
	vnl_vector<double> m_A3;
	vnl_vector<double> m_D1;
	vnl_vector<double> m_D2;
	vnl_vector<double> m_D3;
	vnl_vector<double> m_D4;
	
	//temp container used by some methods
	vnl_matrix<double> m_Container;
	vnl_vector<double> m_Temp;
public:
	/************************* get and set functions *************************/
	void SetNumberOfBins(int aS);
	int GetNumberOfBins();
	void SetRepeat(int aS);
	int GetRepeat();
	void SetSampeSize(int aS);
	int GetSampleSize();
	void SetReference(vnl_vector<double> * aRef);
	vnl_vector<double> * GetReference();
	
	vnl_vector<double>* GetA3();
	vnl_vector<double>* GetD1();
	vnl_vector<double>* GetD2();
	vnl_vector<double>* GetD3();
	vnl_vector<double>* GetD4();
	/************************* operation functions *************************/
	void Initial();
	void Initial(int aNb, int aRepeat, int aSampleSize);
	void ComputeA3();
	void ComputeD1();
	void ComputeD2();
	void ComputeD3();
	void ComputeD4();
	void Compute();
	double Distance(vnl_vector<double> * aPDF1, vnl_vector<double> *aPDF2);
	
		
};

#ifndef ITK_MANUAL_INSTANTIATION
#include "Distribution.hxx"
#endif
	
#endif	
	
