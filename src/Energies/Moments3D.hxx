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

/*
 * 3D version of Moments

 */

#ifndef __Moments3D_hxx
#define __Moments3D_hxx
#include<iostream>

#include "Moments3D.h"

void Moments3D::SetReference(LUMatrix *aRef)
{
    //deep copy
    m_Reference=(*aRef);
    return;
}

LUMatrix * Moments3D::GetReference()
{
    return &m_Reference;
}

LUMatrix* Moments3D::GetMoments()
{
    return &m_Moments;
}

void Moments3D::SetCurrentMeasure(double aM)
{
    m_CurrentMeasure=aM;
}

double Moments3D::GetCurrentMeasure()
{
    return m_CurrentMeasure;
}

void Moments3D::Init()
{
    //set memory
    unsigned int vOrderPlus=GetOrder()+1;
    m_Reference.set_size(vOrderPlus);
    m_Moments.set_size(vOrderPlus);
    m_OnePointMoments.set_size(vOrderPlus);
    m_UpdatedCentralMoments.set_size(vOrderPlus);
    m_UpdatedNormCentralMoments.set_size(vOrderPlus);
    m_UpdatedMoments.set_size(vOrderPlus);
    m_AddPointsMoments.set_size(vOrderPlus);
    m_DeletePointsMoments.set_size(vOrderPlus);
    m_CentralMoments.set_size(vOrderPlus);
    m_NCentralMoments.set_size(vOrderPlus);
    m_Co_Legendre.set_size(vOrderPlus,vOrderPlus);

    m_TempForCentral.set_size(vOrderPlus);
    m_TempForLegendre.set_size(vOrderPlus);
    m_TempForNorm.set_size(vOrderPlus);
    m_TempForSin.set_size(vOrderPlus);
    m_TempForCos.set_size(vOrderPlus);
    m_dxdydz.set_size(vOrderPlus);
    m_dx.set_size(vOrderPlus);
    m_dy.set_size(vOrderPlus);
    m_dz.set_size(vOrderPlus);
    m_NChooseK.set_size(vOrderPlus,vOrderPlus);
    m_GaussianWeight.set_size(vOrderPlus);

    unsigned int vImageSize=GetImageSize();
    m_Co_Chebyshev.set_size(vOrderPlus,vOrderPlus);
    m_Table_Chebyshev.set_size(vOrderPlus,vImageSize);
    m_Table_Legendre.set_size(vOrderPlus,vImageSize);

    //compute Co
    Co_Legendre();
    NChooseK();
    GaussianWeight(200);
}

void Moments3D::PointPosition(unsigned int aX, unsigned int aY, unsigned int aZ, double *aNX, double *aNY, double *aNZ)
{
    //compute point position in the new coordinate
    double cX=GetCenterXDouble();
    double cY=GetCenterYDouble();
    double cZ=GetCenterZDouble();
    double vS=GetScaleDouble();
    *aNX=((double)aX-cX)/vS;
    *aNY=((double)aY-cY)/vS;
    *aNZ=((double)aZ-cZ)/vS;
}

void Moments3D::NChooseK()
{
    unsigned int vOrderPlus=GetOrder()+1;
    //Coefficients
    m_NChooseK.fill(0);
    m_NChooseK(0,0)=1;
    m_NChooseK(1,0)=1;
    m_NChooseK(1,1)=1;
    for (unsigned int i=2; i<vOrderPlus; i++)
    {
        m_NChooseK(i,0)=1;
        for (unsigned int j=1; j<i; j++)
        {
            m_NChooseK(i,j)=m_NChooseK(i-1,j-1)+m_NChooseK(i-1,j);
        }
        m_NChooseK(i,i)=1;
    }

    //std::cout<<m_NChooseK<<std::endl;
}

void Moments3D::Co_Legendre()
{
    //compute legendre co
    //the same with the matlab code
    unsigned int vOrderPlus=GetOrder()+1;
    m_Co_Legendre(0,0)=1;
    for (unsigned int i=1; i<vOrderPlus; i++)
    {
        m_Co_Legendre(i,0)=m_Co_Legendre(i-1,0)*(2*i-1)/i;
        for (unsigned int j=2; j<=i; j++,j++)
        {
            m_Co_Legendre(i,j)=-m_Co_Legendre(i,j-2)*(i-j+2)*(i-j+1)/((2*i-j+1)*j);
        }
    }

    //reverse the coeff
    double temp=0;
    for (unsigned int i=1; i<vOrderPlus; i++)
    {
        for (unsigned int j=0; j<=(i/2); j++)
        {
            temp=m_Co_Legendre(i,j);
            m_Co_Legendre(i,j)=m_Co_Legendre(i,i-j);
            m_Co_Legendre(i,i-j)=temp;
        }
    }
    //std::cout<<m_Co_Legendre<<std::endl;
}

void Moments3D::CopyMoments(LUMatrix *src, LUMatrix *dst)
{
    double *** pS=src->m_Matrix;
    double *** pD=dst->m_Matrix;
    unsigned int vOrderPlus=GetOrder()+1;
    for (unsigned int i=0;i<vOrderPlus;i++)
        for (unsigned int j=0;j<vOrderPlus-i;j++)
            for (unsigned int k=0;k<vOrderPlus-i-j;k++)
            {
                pD[i][j][k]=pS[i][j][k];
            }
}

void Moments3D::NewCenter(double aX, double aY, double aZ, bool AddOrDelete, double *dx, double * dy, double *dz)
{
    //add or delete a point, compute the dx dy
    unsigned int dataS=GetDataSize();
//    double totalx=GetCenterXDouble()*dataS;
//    double totaly=GetCenterYDouble()*dataS;
    if (AddOrDelete)
    {
        *dx=-aX/(dataS+1);
        *dy=-aY/(dataS+1);
        *dz=-aZ/(dataS+1);
    }
    else
    {
        *dx=aX/(dataS-1);
        *dy=aY/(dataS-1);
        *dz=aZ/(dataS-1);
    }
}

double Moments3D::AddOnePoint(unsigned int aX, unsigned int aY, unsigned int aZ)
{
    double nX, nY, nZ, dx, dy, dz;
    PointPosition(aX,aY,aZ,&nX,&nY,&nZ);
    NewCenter(nX,nY,nZ,true,&dx, &dy, &dz);
    TranslateMoments(&m_CentralMoments, &m_UpdatedCentralMoments, dx, dy, dz);
    OnePointMoments(nX+dx, nY+dy, nZ+dz);
    AddOrDeleteMoments(&m_UpdatedCentralMoments, &m_OnePointMoments, true);
    //NormalizeMoments(&m_UpdatedCentralMoments, &m_UpdatedNormCentralMoments);
    ScaleMoments(&m_UpdatedCentralMoments,&m_UpdatedNormCentralMoments);
    LegendreMoments(&m_UpdatedNormCentralMoments,&m_UpdatedMoments);
    m_MomentsMeasure=Distance(&m_UpdatedMoments,&m_Reference, false)-m_CurrentMeasure;
    return m_MomentsMeasure;

}

double Moments3D::DeleteOnePoint(unsigned int aX, unsigned int aY, unsigned int aZ)
{
    double nX, nY, nZ, dx, dy, dz;
    PointPosition(aX,aY,aZ,&nX,&nY,&nZ);
    OnePointMoments(nX, nY, nZ);
    m_TempForCentral=m_CentralMoments;
    //std::cout<<"central moments: "<<m_CentralMoments<<std::endl;
    AddOrDeleteMoments(&m_TempForCentral, &m_OnePointMoments, false);
    //std::cout<<"delete moments: "<<m_TempForCentral<<std::endl;
    NewCenter(nX,nY,nZ,false,&dx, &dy,&dz);
    TranslateMoments(&m_TempForCentral, &m_UpdatedCentralMoments,dx, dy, dz);
    //std::cout<<"update central moments: "<<m_UpdatedCentralMoments<<std::endl;
    //NormalizeMoments(&m_UpdatedCentralMoments, &m_UpdatedNormCentralMoments);
    ScaleMoments(&m_UpdatedCentralMoments,&m_UpdatedNormCentralMoments);
    //std::cout<<"update norm central moments: "<<m_UpdatedNormCentralMoments<<std::endl;
    LegendreMoments(&m_UpdatedNormCentralMoments,&m_UpdatedMoments);
    //std::cout<<"updated Legendre moments: "<<m_UpdatedMoments<<std::endl;
    //std::cout<<"reference moments: "<<m_Reference<<std::endl;
    m_MomentsMeasure=Distance(&m_UpdatedMoments,&m_Reference, false)-m_CurrentMeasure;


    return m_MomentsMeasure;
}

void Moments3D::TranslateMoments(LUMatrix *src, LUMatrix * dst, double dx, double dy, double dz)
{
    //only used for central moments when the center moves
    //x'=x+dx, y'=y+dx


    //first, fill dx^p, dy^q
    unsigned int vOrderPlus=GetOrder()+1;
    m_dx(0)=1;
    m_dy(0)=1;
    m_dz(0)=1;
    for (unsigned int i=1;i<vOrderPlus;i++)
    {
        m_dx(i)=m_dx(i-1)*dx;
        m_dy(i)=m_dy(i-1)*dy;
        m_dz(i)=m_dz(i-1)*dz;
    }
    //second, fill dxdydz
    double temp;
    double ***pDxDyDz=m_dxdydz.m_Matrix;
    for (unsigned int i=0;i<vOrderPlus;i++)
    {
        for (unsigned int j=0;j<vOrderPlus-i;j++)
        {
            temp=m_dx(i)*m_dy(j);
            for (unsigned int k=0;k<vOrderPlus-i-j;k++)
            {
                pDxDyDz[i][j][k]=temp*m_dz(k);
            }
        }
    }
    //third, compute dst
    double *** pSrc=src->m_Matrix;
    double *** pDst=dst->m_Matrix;
    double ** pNChooseK=m_NChooseK.data_array();

    /*
     the method here should be timed and changed to faster algorithm.
     */
    //quite slow in this way, find new algorithm to do so
    clock_t ca, cb;
    ca=clock();
    double temp;
    for (unsigned int i=0;i<vOrderPlus;i++)
        for (unsigned int j=0;j<vOrderPlus-i;j++)
            for (unsigned int w=0;w<vOrderPlus-i-j;w++)
            {
                temp=0;
                for (unsigned int k=0;k<=i;k++)
                    for (unsigned int m=0;m<=j;m++)
                        for (unsigned int ww=0;ww<=w;ww++)
                        {
                            temp+=(pNChooseK[i][k]*pNChooseK[j][m]*pNChooseK[w][ww]*pDxDyDz[i-k][j-m][w-ww]*pSrc[k][m][ww]);
                        }
                pDst[i][j][k]=temp;
            }
    cb=clock();
    //std::cout<<"Translate Moments time: "<<cb-ca<<std::endl;
}

void Moments3D::NormalizeMoments(LUMatrix * src, LUMatrix * dst)
{
    unsigned int vOrderPlus=GetOrder()+1;
    double *** pS=src->m_Matrix;
    double *** pD=dst->m_Matrix;

    double vBase=pS[0][0][0];
    if (vBase<0)
    {
        std::cout<<"Error in NormalizeMoments: the area is negtive."<<std::endl;
    }
    double vRoot=cbrt(vBase);
    //fill m_TempForNorm
    m_TempForNorm(0)=vBase;
    for (unsigned int i=1;i<vOrderPlus;i++)
    {
        m_TempForNorm(i)=m_TempForNorm(i-1)*vRoot;
    }
    //normalize the moments
    for (unsigned int i=0;i<vOrderPlus;i++)
        for (unsigned int j=0;j<vOrderPlus-i;j++)
            for (unsigned int k=0;k<vOrderPlus-i-j;k++)
                pD[i][j][k]=pS[i][j][k]/m_TempForNorm[i+j+k];
}

void Moments3D::ScaleMoments(LUMatrix *src, LUMatrix *dst)
{
    ScaleMoments(src,dst,0.25);
}

void Moments3D::ScaleMoments(LUMatrix *src, LUMatrix *dst, double rate)
{
    //scale the moments to fixed scale, 1/4 of image size by default
    unsigned int vOrderPlus=GetOrder()+1;
    unsigned int vImageSize=GetImageSize();
    double *** pSrc=src->m_Matrix;
    double *** pDst=dst->m_Matrix;
    double area=pSrc[0][0][0];
    double target=rate*(double)vImageSize*(double)vImageSize*(double)vImageSize;
    double root=cbrt(target/area);

    //fill m_TempForNorm
    m_TempForNorm(0)=target/area;
    for (unsigned int i=1;i<vOrderPlus;i++)
    {
        m_TempForNorm(i)=m_TempForNorm(i-1)*root;
    }
    //scale the moments
    for (unsigned int i=0;i<vOrderPlus;i++)
        for (unsigned int j=0;j<vOrderPlus-i;j++)
            for (unsigned int k=0;k<vOrderPlus-i-j;k++)
            {
                pDst[i][j][k]=pSrc[i][j][k]*m_TempForNorm(i+j+k);
            }
}

void Moments3D::RotateMoments(LUMatrix * src, LUMatrix * dst, double aAngle)
{
    //rotate is indenpendent with scale
    double *** pSrc=src->m_Matrix;
    double *** pDst=dst->m_Matrix;
    unsigned int vOrderPlus=GetOrder()+1;
    //compute the sin and cos
    double sint=sin(aAngle);
    double cost=cos(aAngle);
    m_TempForSin(0)=1;
    m_TempForCos(0)=1;
    for (unsigned int i=1;i<vOrderPlus;i++)
    {
        m_TempForSin(i)=m_TempForSin(i-1)*sint;
        m_TempForCos(i)=m_TempForCos(i-1)*cost;
    }
    //rotate in 3D should be handled in different way
    //not ready yet
    /*
     */
}

void Moments3D::LegendreMoments(LUMatrix *GeoMoments, LUMatrix *Legendre)
{
    //compute Legendre Moments from geo moments
    //not ready yet
    unsigned int vOrderPlus=GetOrder()+1;
    double vScale3=GetScaleDouble()*GetScaleDouble()*GetScaleDouble();
    double *** pG=GeoMoments->m_Matrix;
    double *** pL=Legendre->m_Matrix;
    double ** pC=m_Co_Legendre.data_array();
    double *** pT=m_TempForLegendre.m_Matrix;

    //how to compute the moments?
    unsigned int s=0;
    for (unsigned int i=0;i<vOrderPlus;i++)
        for (unsigned int j=0;j<vOrderPlus;j++)
        {
            //pick up the smaller one
            s=i>(vOrderPlus-j-1)?(vOrderPlus-j-1):i;
            for (unsigned int k=0;k<=s;k++)
                pT[i][j]+=(pC[i][k]*pG[k][j]);
        }
    (*Legendre).fill(0);

    double temp;
    for (unsigned int i=0;i<vOrderPlus;i++)
        for (unsigned int j=0;j<vOrderPlus-i;j++)
        {
            s=j;
            for (unsigned int k=0;k<=j;k++)
                pL[i][j]+=(pT[i][k]*pC[j][k]);
            temp=(2*i+1.0)*(2*j+1.0)/8.0/vScale3;
            pL[i][j]*=temp;
        }
    //std::cout<<*GeoMoments<<std::endl;
    //std::cout<<*Legendre<<std::endl;
}

void Moments3D::OnePointMoments(double aX, double aY, double aZ)
{
    unsigned vOrderPlus=GetOrder()+1;
    double *** pM=m_OnePointMoments.m_Matrix;
    pM[0][0][0]=1;
    //fill the first line along x
    for (unsigned int i=1;i<vOrderPlus;i++)
        pM[i][0][0]=pM[i-1][0][0]*aX;
    //fill the first plane
    for (unsigned int i=1;i<vOrderPlus;i++)
        for (unsigned int j=0;j<vOrderPlus-i;j++)
            pM[i][j][0]=pM[i][j-1][0]*aY;
    //fill the all matrix
    for (unsigned int k=1;k<vOrderPlus;k++)
        for (unsigned int i=0;i<vOrderPlus-k;i++)
            for (unsigned int j=0;j<vOrderPlus-k-i;j++)
                pM[i][j][k]=pM[i][j][k-1]*aZ;
}

void Moments3D::AddOrDeleteMoments(LUMatrix *Moments, LUMatrix *OnePoint, bool AddOrDelete)
{
    double *** pM=Moments->m_Matrix;
    double *** pOne=OnePoint->m_Matrix;
    unsigned int vOrderPlus=GetOrder()+1;
    if (AddOrDelete)
    {
        for (unsigned int i=0;i<vOrderPlus;i++)
            for (unsigned int j=0;j<vOrderPlus-i;j++)
                for (unsigned int k=0;k<vOrderPlus-i-j;k++)
                    pM[i][j][k]+=pOne[i][j][k];
    }
    else
    {
        for (unsigned int i=0;i<vOrderPlus;i++)
            for (unsigned int j=0;j<vOrderPlus-i;j++)
                for (unsigned int k=0;k<vOrderPlus-i-j;k++)
                    pM[i][j][k]-=pOne[i][j][k];
    }
}

void Moments3D::GaussianWeight(double aSigma)
{
    unsigned int vOrderPlus=GetOrder()+1;
    double *** pGaussian=m_GaussianWeight.m_Matrix;
    double vCoefficiant=power(1.0/(sqrt(2*ShapePrior::PI)*aSigma),3);
    double temp=2*aSigma*aSigma;
    for (unsigned int i=0; i<vOrderPlus; i++)
    {
        for (unsigned int j=0; j<vOrderPlus-i; j++)
        {
            for (unsigned int k=0; k<vOrderPlus-i-j;k++)
            {
                pGaussian[i][j][k]=vCoefficiant*exp(-(i*i+j*j+k*k)/temp);
            }
        }
    }
}

double Moments3D::Distance(LUMatrix *mom1,LUMatrix *mom2, bool WeightedOrNot)
{
    double *** pM1=mom1->m_Matrix;
    double *** pM2=mom2->m_Matrix;
    double *** pW=m_GaussianWeight.m_Matrix;

    double temp=0;
    double innerTemp=0;
    unsigned int vOrderPlus=GetOrder()+1;
    if (WeightedOrNot)
    {
        for (unsigned int i=0;i<vOrderPlus;i++)
            for (unsigned int j=0;j<vOrderPlus-i;j++)
                for (unsigned int k=0;k<vOrderPlus-i-j;k++)
                {
                    innerTemp=pM1[i][j][k]-pM2[i][j][k];
                    temp+=(innerTemp*innerTemp*pW[i][j][k]);
                }
    }
    else
    {
        for (unsigned int i=0;i<vOrderPlus;i++)
            for (unsigned int j=0;j<vOrderPlus-i;j++)
                for (unsigned int k=0;k<vOrderPlus-i-j;k++)
                {
                    innerTemp=pM1[i][j][k]-pM2[i][j][k];
                    temp+=(innerTemp*innerTemp);
                }
    }

    return temp;
}

void Moments3D::Legendre()
{
    CentralMomentsForC();
    //std::cout<<"Central Moments: "<<m_CentralMoments<<std::endl;
    //NormalizeMoments(&m_CentralMoments,&m_NCentralMoments);
    ScaleMoments(&m_CentralMoments,&m_NCentralMoments);
    //std::cout<<m_NCentralMoments<<std::endl;
    LegendreMoments(&m_NCentralMoments,&m_Moments);
    double t=Distance(&m_Moments, &m_Reference, false);
    SetCurrentMeasure(t);
    //std::cout<<"Moments: "<<m_Moments<<std::endl;
    std::cout<<m_CurrentMeasure<<std::endl;
}

void Moments3D::StepUpdate(vnl_matrix<unsigned int> *aAdd, vnl_matrix<unsigned int> *aDelete)
{
    //step update the geo moments and recompute current moments
    unsigned int AddSize=aAdd->rows();
    unsigned int DeleteSize=aDelete->rows();
    unsigned int vS=GetDataSize();
    double vScale=GetScaleDouble();
    double cX=GetCenterXDouble();
    double cY=GetCenterYDouble();

    //temp var and pointers
    int tX=0, tY=0;
    double dx=0, dy=0;
    unsigned int ** pA=aAdd->data_array();
    unsigned int ** pD=aDelete->data_array();

    if (AddSize==0 && DeleteSize==0)
    {
        std::cout<<"current moments measure: "<<m_CurrentMeasure<<std::endl;
        return;
    }
    if (AddSize==0)
    {
        //only delete
        for (unsigned int i=0;i<DeleteSize;i++)
        {
            tX+=pD[i][0];
            tY+=pD[i][1];
        }

        dx=((double)tX)/(vS-DeleteSize);
        dy=((double)tY)/(vS-DeleteSize);
        PointsMoments(aDelete, false);
        AddOrDeleteMoments(&m_CentralMoments, &m_DeletePointsMoments, false);
        TranslateMoments(&m_CentralMoments,&m_TempForCentral, dx/vScale, dy/vScale);
        m_CentralMoments=m_TempForCentral;
        LegendreMoments(&m_CentralMoments, &m_Moments);
        SetCenterXDouble(cX-dx);
        SetCenterYDouble(cY-dy);
        SetCurrentMeasure(Distance(&m_Moments, &m_Reference, false));
        std::cout<<"current moments measure: "<<m_CurrentMeasure<<std::endl;
        return;
    }
    if (DeleteSize==0)
    {
        //only add
        for (unsigned int i=0;i<AddSize;i++)
        {
            tX+=pA[i][0];
            tY+=pA[i][1];
        }

        dx=-((double)tX)/(vS+AddSize);
        dy=-((double)tY)/(vS+AddSize);
        PointsMoments(aAdd, true);
        AddOrDeleteMoments(&m_CentralMoments, &m_AddPointsMoments, true);
        TranslateMoments(&m_CentralMoments,&m_TempForCentral, dx/vScale, dy/vScale);
        m_CentralMoments=m_TempForCentral;
        LegendreMoments(&m_CentralMoments, &m_Moments);
        SetCenterXDouble(cX-dx);
        SetCenterYDouble(cY-dy);
        SetCurrentMeasure(Distance(&m_Moments, &m_Reference, false));
        std::cout<<"current moments measure: "<<m_CurrentMeasure<<std::endl;
        return;
    }

    //add and delete

    for (unsigned int i=0;i<AddSize;i++)
    {
        tX+=pA[i][0];
        tY+=pA[i][1];
    }
    for (unsigned int i=0;i<DeleteSize;i++)
    {
        tX-=pD[i][0];
        tY-=pD[i][1];
    }
    dx=cX-(cX*vS+(double)tX)/(vS+AddSize-DeleteSize);
    dy=cY-(cY*vS+(double)tY)/(vS+AddSize-DeleteSize);
    PointsMoments(aAdd, true);
    PointsMoments(aDelete, false);
    //std::cout<<m_AddPointsMoments<<std::endl;
    //std::cout<<m_DeletePointsMoments<<std::endl;
    //std::cout<<m_CentralMoments<<std::endl;
    //delete points moments
    AddOrDeleteMoments(&m_CentralMoments, &m_DeletePointsMoments, false);
    //move the rest
    TranslateMoments(&m_CentralMoments,&m_TempForCentral, dx/vScale, dy/vScale);
    m_CentralMoments=m_TempForCentral;
    //move add moments
    TranslateMoments(&m_AddPointsMoments,&m_TempForCentral, dx/vScale, dy/vScale);
    //add the result to central moments
    AddOrDeleteMoments(&m_CentralMoments, &m_TempForCentral, true);
    //std::cout<<m_CentralMoments<<std::endl;
    LegendreMoments(&m_CentralMoments, &m_Moments);
    //update center
    SetCenterXDouble(cX-dx);
    SetCenterYDouble(cY-dy);
    //update current measure
    SetCurrentMeasure(Distance(&m_Moments, &m_Reference, false));
    std::cout<<"current moments measure: "<<m_CurrentMeasure<<std::endl;
    return;
}

void Moments3D::PointsMoments(vnl_matrix<unsigned int> *aPoints, bool AddOrDelete)
{
    unsigned int vOrderPlus=GetOrder()+1;
    unsigned int vS=aPoints->rows();
    double vScale=GetScaleDouble();
    double cX=GetCenterXDouble();
    double cY=GetCenterYDouble();
    double cZ=GetCenterZDouble();

    if (vS==0)
    {
        if (AddOrDelete)
            m_AddPointsMoments.fill(0);
        else
            m_DeletePointsMoments.fill(0);
        return;
    }

    //pointer to the data
    unsigned int ** pData=aPoints->data_array();

    //X, Y, Z
    vnl_vector<double> vX(vS), vY(vS), vZ(vS);
    for (unsigned int i=0;i<vS;i++)
    {
        vX(i)=(pData[i][0]-cX)/vScale;
        vY(i)=(pData[i][1]-cY)/vScale;
        vZ(i)=(pData[i][2]-cZ)/vScale;
    }

    //Xp, yq
    vnl_matrix<double> Xp(vOrderPlus, vS),Yq(vS,vOrderPlus), Zr(vS, vOrderPlus);
    double ** pXp=Xp.data_array();
    double ** pYq=Yq.data_array();
    double ** pZr=Zr.data_array();

    for (unsigned int i=0;i<vS;i++)
    {
        pXp[0][i]=1;
        pYq[i][0]=1;
        pZr[i][0]=1;
    }
    for (unsigned int i=1;i<vOrderPlus;i++)
        for (unsigned int j=0;j<vS;j++)
        {
            pXp[i][j]=pXp[i-1][j]*vX(j);
            pYq[j][i]=pYq[j][i-1]*vY(j);
            pZr[j][i]=pZr[j][i-1]*vZ(j);
        }
    //compute centeral moments
    double *** pM;
    if (AddOrDelete)
    {
        pM=m_AddPointsMoments.m_Matrix;
    }else
    {
        pM=m_DeletePointsMoments.m_Matrix;
    }
    double temp;
    for (unsigned int i=0;i<vOrderPlus;i++)
        for (unsigned int j=0;j<vOrderPlus-i;j++)
            for (unsigned int k=0;k<vOrderPlus-i-j;k++)
            {
                temp=0;
                for (unsigned w=0;w<vS;w++)
                {
                    temp+=pXp[i][w]*pYq[w][j]*pZr[w][k];
                }
                pM[i][j][k]=temp;
            }

}
void Moments3D::CentralMomentsForC()
{
    //size, order, scale, center
    unsigned int vS=GetDataSize();
    unsigned int vOrderPlus=GetOrder()+1;
    double vScale=GetScaleDouble();
    double cX=GetCenterXDouble();
    double cY=GetCenterYDouble();
    double cZ=GetCenterZDouble();

    //pointer to the data
    unsigned int ** pData=GetInputData()->data_array();

    //X, Y, Z
    vnl_vector<double> vX(vS), vY(vS), vZ(vS);
    for (unsigned int i=0;i<vS;i++)
    {
        vX(i)=(pData[i][0]-cX)/vScale;
        vY(i)=(pData[i][1]-cY)/vScale;
        vZ(i)=(pData[i][2]-cZ)/vScale;
    }

    //Xp, yq
    vnl_matrix<double> Xp(vOrderPlus, vS),Yq(vS,vOrderPlus), Zr(vS, vOrderPlus);
    double ** pXp=Xp.data_array();
    double ** pYq=Yq.data_array();
    double ** pZr=Zr.data_array();

    for (unsigned int i=0;i<vS;i++)
    {
        pXp[0][i]=1;
        pYq[i][0]=1;
        pZr[i][0]=1;
    }
    for (unsigned int i=1;i<vOrderPlus;i++)
        for (unsigned int j=0;j<vS;j++)
        {
            pXp[i][j]=pXp[i-1][j]*vX(j);
            pYq[j][i]=pYq[j][i-1]*vY(j);
            pZr[j][i]=pZr[j][i-1]*vZ(j);
        }
    //compute centeral moments
    double *** pCentral=m_CentralMoments.m_Matrix;
    double temp;
    for (unsigned int i=0;i<vOrderPlus;i++)
        for (unsigned int j=0;j<vOrderPlus-i;j++)
            for (unsigned int k=0;k<vOrderPlus-i-j;k++)
            {
                temp=0;
                for (unsigned w=0;w<vS;w++)
                {
                    temp+=pXp[i][w]*pYq[w][j]*pZr[w][k];
                }
                pM[i][j][k]=temp;
            }

    //std::cout<<m_CentralMoments<<std::endl;
}

//void Moments3D::ComputeTableChebyshev(unsigned int aSize, unsigned int aOrder)
//{
//
//    for (unsigned int j=0;j<aSize;j++)
//    {
//        m_Table_Chebyshev(0,j)=1;
//        m_Table_Chebyshev(1,j)=2.0/(double)aSize*j+(1.0-(int)aSize)/(double)aSize;
//
//    }
//
//    double vNN=(double)aSize;
//    double k, s;
//    for(unsigned int i=2;i<=aOrder;i++)
//        for (unsigned int j=0;j<aSize;j++)
//        {
//            k=(double)i;
//            s=(double)j;
//            m_Table_Chebyshev(i,j)=((2*k-1)/k)*((2*s+1-vNN)/vNN)*m_Table_Chebyshev(i-1,j)-(k-1)/k*(1-((k-1)/vNN)*((k-1)/vNN))*m_Table_Chebyshev(i-2,j);
//        }
//}
//void Moments3D::ChebyshevMoments()
//{
//    clock_t ca,cb;
//    ca=clock();
//    unsigned int **pData=GetInputData()->data_array();
//    double **pMoments=m_Moments.data_array();
//    double **pTable=m_Table_Chebyshev.data_array();
//    unsigned int vSize=GetDataSize();
//    unsigned int vOrder=GetOrder();
//    for (unsigned int i=0;i<vSize;i++)
//    {
//        for(unsigned int k=0;k<=vOrder;k++)
//            for (unsigned int m=0;m<=vOrder;m++)
//            {
//                pMoments[k][m]+=(pTable[k][pData[i][0]]*pTable[m][pData[i][1]]);
//            }
//    }
//    cb=clock();
//    std::cout<<"compute Chebyshev Moments: "<<cb-ca<<std::endl;
//    ca=clock();
//    for(unsigned int k=0;k<=vOrder;k++)
//            for (unsigned int m=0;m<=vOrder;m++)
//            {
//                pMoments[k][m]+=(pTable[k][pData[1][0]]*pTable[m][pData[1][1]]);
//            }
//    cb=clock();
//    std::cout<<"one point moments: "<<cb-ca<<std::endl;
//}


#endif
