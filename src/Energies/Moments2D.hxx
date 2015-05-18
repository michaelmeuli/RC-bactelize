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
 * in this file, several kinds of moments are given as well as the updating method
 * which is used for moments updating during active contour evolution. Please remember
 * that how to compute the moments does not impact the performace while how to update
 * does. However, different moments have different properties which make the evolution
 * different.
 *
 * for input data, scale and translation have been handled at beginning. So, all
 * moments do not have to consider image resolution problem while object's size still
 * should be handled by moments.
 *
 * Be very careful with different moments with the same name. So fa as I know, there
 * are five different kinds of moments named as Chebyshev moments or Tchebichef moments.
 * Even for the same moments, there are a lot of methods to compute them. For this file,
 * due to update, moments are computed directly from geometric moments.

 */

#ifndef __Moments2D_hxx
#define __Moments2D_hxx
#include<iostream>

#include "Moments2D.h"
#include "CellListedHashMap.h"
namespace itk{

void Moments2D::SetReference(vnl_matrix<double> *aRef)
{
    //deep copy
    m_Reference=(*aRef);
    return;
}

vnl_matrix<double> * Moments2D::GetReference()
{
    return &m_Reference;
}

void Moments2D::ComputeReference()
{
    CentralMomentsForC();
    vnl_matrix_fixed<double,3,3> rotation_matrix = RotationMatrixFromMoments(&m_CentralMoments);
    SetPriorAxesVector(&rotation_matrix);
    //store the reference rotation matrix  
    //std::cout<<"the reference axis is: "<<rotation_matrix(0,0)<<" "<<rotation_matrix(1,0)<<" "<<rotation_matrix(0,1)<<" "<<rotation_matrix(1,1)<<" "<<std::endl;    ScaleMoments(&m_CentralMoments,&m_NCentralMoments);
    CopyMoments(&m_NCentralMoments, &m_ReferenceNCentralMoments);
    LegendreMoments(&m_NCentralMoments,&m_Reference);  
}

vnl_matrix<double>* Moments2D::GetMoments()
{
    return &m_Moments;
}

void Moments2D::SetCurrentMeasure(double aM)
{
    m_CurrentMeasure=aM;
}

double Moments2D::GetCurrentMeasure()
{
    return m_CurrentMeasure;
}

void Moments2D::Init()
{
    /*******  set memory   *******/
    unsigned int vOrderPlus=GetOrder()+1;
    
    m_Reference.set_size(vOrderPlus,vOrderPlus);
    m_ReferenceNCentralMoments.set_size(vOrderPlus,vOrderPlus);
    m_Moments.set_size(vOrderPlus,vOrderPlus);
    m_OnePointMoments.set_size(vOrderPlus,vOrderPlus);
    m_UpdatedCentralMoments.set_size(vOrderPlus,vOrderPlus);
    m_UpdatedNormCentralMoments.set_size(vOrderPlus,vOrderPlus);
    m_UpdatedMoments.set_size(vOrderPlus,vOrderPlus);
    m_AddPointsMoments.set_size(vOrderPlus,vOrderPlus);
    m_DeletePointsMoments.set_size(vOrderPlus,vOrderPlus);
    m_CentralMoments.set_size(vOrderPlus,vOrderPlus);
    m_NCentralMoments.set_size(vOrderPlus,vOrderPlus);
    m_RCentralMoments.set_size(vOrderPlus,vOrderPlus);
    
    //coefficient table for continuous moments
    m_Co_Legendre.set_size(vOrderPlus,vOrderPlus);
    m_Co_Chebyshev.set_size(vOrderPlus,vOrderPlus);
    m_NChooseK.set_size(vOrderPlus,vOrderPlus);
    m_GaussianWeight.set_size(vOrderPlus,vOrderPlus);

    //temp var
    m_TempForCentral.set_size(vOrderPlus,vOrderPlus);
    m_TempForRotation.set_size(vOrderPlus,vOrderPlus);
    m_TempForLegendre.set_size(vOrderPlus,vOrderPlus);
    m_TempForNorm.set_size(vOrderPlus);
    m_TempForSin.set_size(vOrderPlus);
    m_TempForCos.set_size(vOrderPlus);
    m_dxdy.set_size(vOrderPlus,vOrderPlus);
    m_dx.set_size(vOrderPlus);
    m_dy.set_size(vOrderPlus);
    
    //coefficient table for discrete moments
    unsigned int vImageSize=GetImageSize();
    m_Table_Chebyshev.set_size(vOrderPlus,vImageSize);
    m_Table_Legendre.set_size(vOrderPlus,vImageSize);
    
    //for angle integration
    m_Integrate_Angle.set_size(vOrderPlus, vOrderPlus);
    
    /******* compute Coefficient *******/
    //only add what is needed here
    Co_Legendre();
    NChooseK();
    Co_Integrate_Angle();
    GaussianWeight(200);
}

void Moments2D::PointPosition(unsigned int aX, unsigned int aY, double *aNX, double *aNY)
{
    //compute point position in the new coordinate
    double cX=GetCenterXDouble();
    double cY=GetCenterYDouble();
    double vS=GetScaleDouble();
    *aNX=((double)aX-cX)/vS;
    *aNY=((double)aY-cY)/vS;
}

void Moments2D::NChooseK()
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

void Moments2D::Co_Legendre()
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

void Moments2D::Co_Integrate_Angle()
{
    //compute the integrate angle table

    //integrate(cos(t)^k*sin(t)^m,t,2,2*pi)
    //2beta((k+1)/2,(m+1)/2) iff k and m are even, zero else

    double **pData=m_Integrate_Angle.data_array();
    m_Integrate_Angle.fill(0);
    
    unsigned int vOrderPlus=GetOrder()+1;
    unsigned int vMaxIndex=vOrderPlus+vOrderPlus;
    //temp container for log-gamma function value
    vnl_vector<double> vLogGammaValue(vMaxIndex);

    for (unsigned int i=0;i<vMaxIndex;i++)
    {
        vLogGammaValue(i)=vnl_log_gamma((i+1)/2.0);
    }
    //integrate table
    for (unsigned int i=0;i<vOrderPlus;i++,i++)
        for (unsigned int j=i;j<vOrderPlus;j++,j++)
        {
            pData[i][j]=2*vcl_exp(vLogGammaValue(i)+vLogGammaValue(j)-vLogGammaValue(i+j+1));
            pData[j][i]=pData[i][j];//the beta function is symmetric
        }
}

void Moments2D::CopyMoments(vnl_matrix<double> *src, vnl_matrix<double> *dst)
{
    double ** pS=src->data_array();
    double ** pD=dst->data_array();
    unsigned int vOrderPlus=GetOrder()+1;
    for (unsigned int i=0;i<vOrderPlus;i++)
        for (unsigned int j=0;j<vOrderPlus-i;j++)
        {
            pD[i][j]=pS[i][j];
        }
}

void Moments2D::NewCenter(double aX, double aY, bool AddOrDelete, double *dx, double * dy)
{
    //add or delete a point, compute the dx dy
    unsigned int dataS=GetDataSize();
//    double totalx=GetCenterXDouble()*dataS;
//    double totaly=GetCenterYDouble()*dataS;
    if (AddOrDelete)
    {
        *dx=-aX/(dataS+1);
        *dy=-aY/(dataS+1);
    }
    else
    {
        *dx=aX/(dataS-1);
        *dy=aY/(dataS-1);
    }
}

double Moments2D::AddOnePoint(unsigned int aX, unsigned int aY)
{
    bool debug = true;
    if (debug)
    {
    ///* back up the old one
    double nX, nY, dx, dy;
    PointPosition(aX,aY,&nX,&nY);
    NewCenter(nX,nY,true,&dx, &dy);
    TranslateMoments(&m_CentralMoments, &m_UpdatedCentralMoments, dx, dy);
    OnePointMoments(nX+dx, nY+dy);
    AddOrDeleteMoments(&m_UpdatedCentralMoments, &m_OnePointMoments, true);
    //NormalizeMoments(&m_UpdatedCentralMoments, &m_UpdatedNormCentralMoments);
    double angle = GetAngleA();
    RotateMoments(&m_UpdatedCentralMoments,&m_OnePointMoments,angle);
    ScaleMoments(&m_OnePointMoments,&m_UpdatedNormCentralMoments);
    LegendreMoments(&m_UpdatedNormCentralMoments,&m_UpdatedMoments);
    m_MomentsMeasure=Distance(&m_UpdatedMoments,&m_Reference, false)-m_CurrentMeasure;
    return m_MomentsMeasure;
    //*/
    }
    
    else
    {
    ///* //this method is trying to avoid rotation for each point. 
    double nX, nY, rX, rY;
    PointPosition(aX,aY,&nX,&nY);
    double vAngle = GetAngleA();
    rX = cos(vAngle)*nX - sin(vAngle)*nY;
    rY = sin(vAngle)*nX + cos(vAngle)*nY;
    
    rX/=m_MomentsScaleRatio;
    rY/=m_MomentsScaleRatio;
    
    
    OnePointMoments(rX, rY);
    AddOrDeleteMoments(&m_RCentralMoments, &m_OnePointMoments, true);
    LegendreMoments(&m_OnePointMoments,&m_UpdatedMoments);
    m_MomentsMeasure=Distance(&m_UpdatedMoments,&m_Reference, false)-m_CurrentMeasure;
    return m_MomentsMeasure;
    //*/
    }
    
}

double Moments2D::DeleteOnePoint(unsigned int aX, unsigned int aY)
{
    bool debug = true;
    if (debug)
    {
    ///*
    double nX, nY, dx, dy;
    PointPosition(aX,aY,&nX,&nY);
    OnePointMoments(nX, nY);
    m_TempForCentral=m_CentralMoments;
    //std::cout<<"central moments: "<<m_CentralMoments<<std::endl;
    AddOrDeleteMoments(&m_TempForCentral, &m_OnePointMoments, false);
    //std::cout<<"delete moments: "<<m_TempForCentral<<std::endl;
    NewCenter(nX,nY,false,&dx, &dy);
    TranslateMoments(&m_TempForCentral, &m_UpdatedCentralMoments,dx, dy);
    //std::cout<<"update central moments: "<<m_UpdatedCentralMoments<<std::endl;
    //NormalizeMoments(&m_UpdatedCentralMoments, &m_UpdatedNormCentralMoments);
    double angle = GetAngleA();
    RotateMoments(&m_UpdatedCentralMoments,&m_OnePointMoments,angle);
    ScaleMoments(&m_OnePointMoments,&m_UpdatedNormCentralMoments);
    //std::cout<<"update norm central moments: "<<m_UpdatedNormCentralMoments<<std::endl;
    LegendreMoments(&m_UpdatedNormCentralMoments,&m_UpdatedMoments);
    //std::cout<<"updated Legendre moments: "<<m_UpdatedMoments<<std::endl;
    //std::cout<<"reference moments: "<<m_Reference<<std::endl;
    m_MomentsMeasure=Distance(&m_UpdatedMoments,&m_Reference, false)-m_CurrentMeasure;
    return m_MomentsMeasure;
    //*/
    }
    else
    {
    ///* //this method is trying to avoid rotation for each point. 
    double nX, nY, rX, rY;
    PointPosition(aX,aY,&nX,&nY);
    double vAngle = GetAngleA();
    rX = cos(vAngle)*nX - sin(vAngle)*nY;
    rY = sin(vAngle)*nX + cos(vAngle)*nY;
    
    rX/=m_MomentsScaleRatio;
    rY/=m_MomentsScaleRatio;
    
    
    OnePointMoments(rX, rY);
    AddOrDeleteMoments(&m_RCentralMoments, &m_OnePointMoments, false);
    LegendreMoments(&m_OnePointMoments,&m_UpdatedMoments);
    m_MomentsMeasure=Distance(&m_UpdatedMoments,&m_Reference, false)-m_CurrentMeasure;
    
    return m_MomentsMeasure;
    //*/
    }
}

void Moments2D::TranslateMoments(vnl_matrix<double> *src, vnl_matrix<double> * dst, double dx, double dy)
{
    //only used for central moments when the center moves
    //x'=x+dx, y'=y+dx


    //first, fill dx^p, dy^q
    unsigned int vOrderPlus=GetOrder()+1;
    m_dx(0)=1;
    m_dy(0)=1;
    for (unsigned int i=1;i<vOrderPlus;i++)
    {
        m_dx(i)=m_dx(i-1)*dx;
        m_dy(i)=m_dy(i-1)*dy;
    }
    //second, fill dxdy
    double **pDxDy=m_dxdy.data_array();
    for (unsigned int i=0;i<vOrderPlus;i++)
        for (unsigned int j=0;j<vOrderPlus-i;j++)
        {
            pDxDy[i][j]=m_dx(i)*m_dy(j);
        }
    //third, compute dst
    double ** pSrc=src->data_array();
    double ** pDst=dst->data_array();
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
        {
            temp=0;
            for (unsigned int k=0;k<=i;k++)
                for (unsigned int m=0;m<=j;m++)
                {
                    temp+=(pNChooseK[i][k]*pNChooseK[j][m]*pDxDy[i-k][j-m]*pSrc[k][m]);
                }
            pDst[i][j]=temp;
        }
    cb=clock();
    //std::cout<<"Translate Moments time: "<<cb-ca<<std::endl;
}

void Moments2D::NormalizeMoments(vnl_matrix<double> * src, vnl_matrix<double> * dst)
{
    unsigned int vOrderPlus=GetOrder()+1;
    double ** pS=src->data_array();
    double ** pD=dst->data_array();

    double vBase=pS[0][0];
    if (vBase<0)
    {
        std::cout<<"Error in NormalizeMoments: the area is negtive."<<std::endl;
    }
    double vRoot=sqrt(vBase);
    //fill m_TempForNorm
    m_TempForNorm(0)=vBase;
    for (unsigned int i=1;i<vOrderPlus;i++)
    {
        m_TempForNorm(i)=m_TempForNorm(i-1)*vRoot;
    }
    //normalize the moments
    for (unsigned int i=0;i<vOrderPlus;i++)
        for (unsigned int j=0;j<vOrderPlus-i;j++)
            pD[i][j]=pS[i][j]/m_TempForNorm[i+j];
}

void Moments2D::ScaleMoments(vnl_matrix<double> *src, vnl_matrix<double> *dst)
{
    ScaleMoments(src,dst,0.5);
}

void Moments2D::ScaleMoments(vnl_matrix<double> *src, vnl_matrix<double> *dst, double rate)
{
    //scale the moments to fixed scale, 1/4 of image size by default
    unsigned int vOrderPlus=GetOrder()+1;
    unsigned int vImageSize=GetImageSize();
    double ** pSrc=src->data_array();
    double ** pDst=dst->data_array();
    double area=pSrc[0][0];
    double target=rate*(double)vImageSize*(double)vImageSize;
    double root=sqrt(target/area);

    //fill m_TempForNorm
    m_TempForNorm(0)=target/area;
    for (unsigned int i=1;i<vOrderPlus;i++)
    {
        m_TempForNorm(i)=m_TempForNorm(i-1)*root;
    }
    //scale the moments
    for (unsigned int i=0;i<vOrderPlus;i++)
        for (unsigned int j=0;j<vOrderPlus-i;j++)
        {
            pDst[i][j]=pSrc[i][j]*m_TempForNorm(i+j);
        }
}

void Moments2D::RotateMoments(vnl_matrix<double> * src, vnl_matrix<double> * dst, double aAngle)
{
    //rotate is indenpendent with scale
    double ** pSrc=src->data_array();
    double ** pDst=dst->data_array();
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
    double temp=0;
    clock_t ca, cb;
    ca=clock();
    for (unsigned int i=0;i<vOrderPlus;i++)
        for (unsigned int j=0;j<vOrderPlus-i;j++)
        {
            temp=0;
            for (unsigned int k=0;k<=i;k++)
                for (unsigned int m=0;m<=j;m++)
                {
                    if ((j-m)%2==0)
                        temp+=(pSrc[j+k-m][i+m-k]*m_NChooseK(i,k)*m_NChooseK(j,m)*m_TempForSin(i+j-k-m)*m_TempForCos(k+m));
                    else
                        temp-=(pSrc[j+k-m][i+m-k]*m_NChooseK(i,k)*m_NChooseK(j,m)*m_TempForSin(i+j-k-m)*m_TempForCos(k+m));
                }
            pDst[i][j]=temp;
        }
    cb=clock();
    //std::cout<<"Rotate Moments time: "<<cb-ca<<std::endl;
}

void Moments2D::IntegrateAngle(vnl_matrix<double> *src, vnl_matrix<double> *dst)
{
    //compute the integration of angle to make rotation invariant
    double ** pS=src->data_array();
    double ** pD=dst->data_array();
    double ** pTableAngle=m_Integrate_Angle.data_array();

    dst->fill(0);
    unsigned int vOrderPlus=GetOrder()+1;

    clock_t ca, cb;
    ca=clock();
    for (unsigned int i=0;i<vOrderPlus;i++)
        for (unsigned int j=0;j<vOrderPlus-i;j++)
        {
            for (unsigned int k=0;k<=i;k++)
                for (unsigned int m=0;m<=j;m++)
                {
                    if ((i-k+m)%2==0 && (j-m+k)%2==0)
                    {
                        if ((i-k)%2==0)
                            pD[i][j]+=pTableAngle[i-k+m][j-m+k]*pS[m+k][i+j-k-m];
                        else
                            pD[i][j]-=pTableAngle[i-k+m][j-m+k]*pS[m+k][i+j-k-m];
                    }
                }
        }
    cb=clock();
    std::cout<<"Integrate Angle Cost Time: "<<cb-ca<<std::endl;
    //notice: the table of angle integration is constant and reduces rapidly when order goes larger.
    //therefore, a proper approximation could be used.
}
void Moments2D::LegendreMoments(vnl_matrix<double> *GeoMoments, vnl_matrix<double> *Legendre)
{
    //compute Legendre Moments from geo moments
    unsigned int vOrderPlus=GetOrder()+1;
    double vScale2=GetScaleDouble()*GetScaleDouble();
    double ** pG=GeoMoments->data_array();
    double ** pL=Legendre->data_array();
    double ** pC=m_Co_Legendre.data_array();
    double ** pT=m_TempForLegendre.data_array();

    m_TempForLegendre.fill(0);
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
            temp=(2*i+1.0)*(2*j+1.0)/4.0/vScale2;
            pL[i][j]*=temp;
        }
    //std::cout<<*GeoMoments<<std::endl;
    //std::cout<<*Legendre<<std::endl;
}

void Moments2D::OnePointMoments(double aX, double aY)
{
    unsigned vOrderPlus=GetOrder()+1;
    double ** pM=m_OnePointMoments.data_array();
    pM[0][0]=1;
    for (unsigned int i=1;i<vOrderPlus;i++)
        pM[0][i]=pM[0][i-1]*aY;
    for (unsigned int i=1;i<vOrderPlus;i++)
        for (unsigned int j=0;j<vOrderPlus-i;j++)
            pM[i][j]=pM[i-1][j]*aX;
}

void Moments2D::AddOrDeleteMoments(vnl_matrix<double> *Moments, vnl_matrix<double> *OnePoint, bool AddOrDelete)
{
    double ** pM=Moments->data_array();
    double ** pOne=OnePoint->data_array();
    unsigned int vOrderPlus=GetOrder()+1;
    if (AddOrDelete)
    {
        for (unsigned int i=0;i<vOrderPlus;i++)
            for (unsigned int j=0;j<vOrderPlus-i;j++)
                pM[i][j]+=pOne[i][j];
    }
    else
    {
        for (unsigned int i=0;i<vOrderPlus;i++)
            for (unsigned int j=0;j<vOrderPlus-i;j++)
                pM[i][j]-=pOne[i][j];
    }
}

void Moments2D::GaussianWeight(double aSigma)
{
    unsigned int vOrderPlus=GetOrder()+1;
    double ** pGaussian=m_GaussianWeight.data_array();
    double vCoefficiant=1.0/(2*ShapePrior::PI*aSigma*aSigma);
    for (unsigned int i=0; i<vOrderPlus; i++)
    {
        for (unsigned int j=0; j<vOrderPlus-i; j++)
        {
            pGaussian[i][j]=vCoefficiant*exp(-(i*i+j*j)/(2*aSigma*aSigma));
        }
    }
}

double Moments2D::Distance(vnl_matrix<double> *mom1,vnl_matrix<double> *mom2, bool WeightedOrNot)
{
    double ** pM1=mom1->data_array();
    double ** pM2=mom2->data_array();
    double ** pW=m_GaussianWeight.data_array();

    double temp=0;
    double innerTemp=0;
    unsigned int vOrderPlus=GetOrder()+1;
    if (WeightedOrNot)
    {
        for (unsigned int i=0;i<vOrderPlus;i++)
            for (unsigned int j=0;j<vOrderPlus-i;j++)
            {
                innerTemp=pM1[i][j]-pM2[i][j];
                temp+=(innerTemp*innerTemp*pW[i][j]);
            }
    }
    else
    {
        for (unsigned int i=0;i<vOrderPlus;i++)
            for (unsigned int j=0;j<vOrderPlus-i;j++)
            {
                innerTemp=pM1[i][j]-pM2[i][j];
                temp+=(innerTemp*innerTemp);
            }
    }

    return sqrt(temp);
}

void Moments2D::Legendre()
{
    bool debug = false;
    if (debug)
    {
    //translation
    CentralMomentsForC();
    
    //scale
    ScaleMoments(&m_CentralMoments,&m_NCentralMoments);
    
    
    //orthogonal moments
    LegendreMoments(&m_NCentralMoments,&m_Moments);
    double t=Distance(&m_Moments, &m_Reference, false);
    SetCurrentMeasure(t);
    //std::cout<<"Moments: "<<m_Moments<<std::endl;
    std::cout<<m_CurrentMeasure<<std::endl;
    }
    else
    {
    //translation
    CentralMomentsForC();
    
    
    //rotation
    m_TempRotationMatrix = RotationMatrixFromMoments(&m_CentralMoments);
    double vAngle = TheBestAngle(&m_CentralMoments, &m_TempRotationMatrix, &m_Reference, GetPriorAxesVector());
    //the negative sign is from the different coord system
    SetAngleA(-vAngle);   
    std::cout<<"current angle"<<-vAngle<<std::endl;
    RotateMoments(&m_CentralMoments, &m_RCentralMoments, vAngle);
    
    //scale
    ScaleMoments(&m_RCentralMoments,&m_NCentralMoments);    
    //orthogonal moments
    LegendreMoments(&m_NCentralMoments,&m_Moments);
    double t=Distance(&m_Moments, &m_Reference, false);
    SetCurrentMeasure(t);
    //std::cout<<"Moments: "<<m_Moments<<std::endl;
    std::cout<<m_CurrentMeasure<<std::endl;
    }
}

void Moments2D::StepUpdate(vnl_matrix<unsigned int> *aAdd, vnl_matrix<unsigned int> *aDelete)
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

void Moments2D::PointsMoments(vnl_matrix<unsigned int> *aPoints, bool AddOrDelete)
{
    unsigned int vOrderPlus=GetOrder()+1;
    unsigned int vS=aPoints->rows();
    double vScale=GetScaleDouble();
    double cX=GetCenterXDouble();
    double cY=GetCenterYDouble();

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

    //X, Y
    vnl_vector<double> vX(vS), vY(vS);
    for (unsigned int i=0;i<vS;i++)
    {
        vX(i)=(pData[i][0]-cX)/vScale;
        vY(i)=(pData[i][1]-cY)/vScale;
    }

    //Xp, yq
    vnl_matrix<double> Xp(vOrderPlus, vS),Yq(vS,vOrderPlus);
    double ** pXp=Xp.data_array();
    double ** pYq=Yq.data_array();

    for (unsigned int i=0;i<vS;i++)
    {
        pXp[0][i]=1;
        pYq[i][0]=1;
    }
    for (unsigned int i=1;i<vOrderPlus;i++)
        for (unsigned int j=0;j<vS;j++)
        {
            pXp[i][j]=pXp[i-1][j]*vX(j);
            pYq[j][i]=pYq[j][i-1]*vY(j);
        }
    //compute centeral moments
    double ** pM;
    if (AddOrDelete)
    {
        pM=m_AddPointsMoments.data_array();
    }else
    {
        pM=m_DeletePointsMoments.data_array();
    }
    double temp;
    for (unsigned int i=0;i<vOrderPlus;i++)
        for (unsigned int j=0;j<vOrderPlus-i;j++)
        {
            temp=0;
            for (unsigned k=0;k<vS;k++)
            {
                temp+=pXp[i][k]*pYq[k][j];
            }
            pM[i][j]=temp;
        }

}
void Moments2D::CentralMomentsForC()
{
    //size, order, scale, center
    unsigned int vS=GetDataSize();
    unsigned int vOrderPlus=GetOrder()+1;
    double vScale=GetScaleDouble();
    double cX=GetCenterXDouble();
    double cY=GetCenterYDouble();

    //pointer to the data
    unsigned int ** pData=GetInputData()->data_array();

    //X, Y
    vnl_vector<double> vX(vS), vY(vS);
    for (unsigned int i=0;i<vS;i++)
    {
        vX(i)=(pData[i][0]-cX)/vScale;
        vY(i)=(pData[i][1]-cY)/vScale;
    }

    //Xp, yq
    vnl_matrix<double> Xp(vOrderPlus, vS),Yq(vS,vOrderPlus);
    double ** pXp=Xp.data_array();
    double ** pYq=Yq.data_array();
    
    for (unsigned int i=0;i<vS;i++)
    {
        pXp[0][i]=1;
        pYq[i][0]=1;
    }
    for (unsigned int i=1;i<vOrderPlus;i++)
        for (unsigned int j=0;j<vS;j++)
        {
            pXp[i][j]=pXp[i-1][j]*vX(j);
            pYq[j][i]=pYq[j][i-1]*vY(j);
        }
    //compute centeral moments
    double ** pCentral=m_CentralMoments.data_array();
    double temp;
    for (unsigned int i=0;i<vOrderPlus;i++)
        for (unsigned int j=0;j<vOrderPlus-i;j++)
        {
            temp=0;
            for (unsigned k=0;k<vS;k++)
            {
                temp+=pXp[i][k]*pYq[k][j];
            }
            pCentral[i][j]=temp;
        }

    //std::cout<<m_CentralMoments<<std::endl;
}

void Moments2D::Table_Chebyshev(unsigned int aSize, unsigned int aOrder)
{
    clock_t ca,cb;
    ca=clock();
    for (unsigned int j=0;j<aSize;j++)
    {
        m_Table_Chebyshev(0,j)=1;
        m_Table_Chebyshev(1,j)=2.0/(double)aSize*j+(1.0-(int)aSize)/(double)aSize;

    }

    double vNN=(double)aSize;
    double k, s;
    for(unsigned int i=2;i<=aOrder;i++)
        for (unsigned int j=0;j<aSize;j++)
        {
            k=(double)i;
            s=(double)j;
            m_Table_Chebyshev(i,j)=((2*k-1)/k)*((2*s+1-vNN)/vNN)*m_Table_Chebyshev(i-1,j)-(k-1)/k*(1-((k-1)/vNN)*((k-1)/vNN))*m_Table_Chebyshev(i-2,j);
        }
    cb=clock();
    std::cout<<"Table Construction time: "<<cb-ca<<std::endl;
}

void Moments2D::ChebyshevMoments(vnl_matrix<unsigned int> * aCoord)
{
    clock_t ca,cb;
    ca=clock();
    unsigned int **pData=aCoord->data_array();
    double **pMoments=m_Moments.data_array();
    double **pTable=m_Table_Chebyshev.data_array();
    unsigned int vSize=GetDataSize();
    unsigned int vOrder=GetOrder();
    for (unsigned int i=0;i<vSize;i++)
    {
        for(unsigned int k=0;k<=vOrder;k++)
            for (unsigned int m=0;m<=vOrder;m++)
            {
                pMoments[k][m]+=(pTable[k][pData[i][0]]*pTable[m][pData[i][1]]);
            }
    }
    cb=clock();
    std::cout<<"compute Chebyshev Moments: "<<cb-ca<<std::endl;
    ca=clock();
    for(unsigned int k=0;k<=vOrder;k++)
            for (unsigned int m=0;m<=vOrder;m++)
            {
                pMoments[k][m]+=(pTable[k][pData[1][0]]*pTable[m][pData[1][1]]);
            }
    cb=clock();
    std::cout<<"one point moments: "<<cb-ca<<std::endl;
}

double Moments2D::TheBestAngle(vnl_matrix<double> *aMom1, vnl_matrix_fixed<double,3,3> *aVectors1, vnl_matrix<double> *aMom2 , vnl_matrix_fixed<double,3,3> *aVectors2)
{
    //attention: the two moments must have been scaled to central moments
    
    //compute the best angle which align the moments

    vnl_vector_fixed<double,2> vect1AxisLong, vect2AxisLong;
    for (unsigned int i = 0; i<2;i++)
    {
        vect1AxisLong[i]=(*aVectors1)[i][0];
        vect2AxisLong[i]=(*aVectors2)[i][0];
    }
    
    std::cout<<"rotation vectors: "<<vect1AxisLong[0]<<" "<<vect1AxisLong[1]<<" "<<vect2AxisLong[0]<<" "<<vect2AxisLong[1]<<std::endl;
    //now compute the best rotation angle
    double sinTheta = (vect2AxisLong[1]*vect1AxisLong[0] - vect1AxisLong[1]*vect2AxisLong[0])/(vect1AxisLong[0]*vect1AxisLong[0]+vect1AxisLong[1]*vect1AxisLong[1]);
    
    double Theta = asin(sinTheta);
    RotateMoments(aMom1,&m_TempForRotation,Theta);
    double dist1 = Distance(aMom2,&m_TempForRotation,false);
    RotateMoments(aMom1,&m_TempForRotation,Theta+ShapePrior::PI);
    double dist2 = Distance(aMom2,&m_TempForRotation,false);
    
    if (dist1<=dist2)
    {
        return Theta;
    }else
    {
        return Theta+ShapePrior::PI;
    }
    
}

vnl_matrix_fixed<double,3,3> Moments2D::RotationMatrixFromMoments(vnl_matrix<double> * CentralMoments)
{
    //compute the rotation matrix from low order moments
    double ** p = CentralMoments->data_array();
    double xx = p[2][0];
    double yy = p[0][2];
    double xy = p[1][1];
    
    vnl_matrix<double> M(2,2);
    M[0][0]=xx;
    M[0][1]=xy;
    M[1][0]=xy;
    M[1][1]=yy;
    
    //vnl_svd_economy<double> mm(M);
    vnl_svd_economy<double> mm(M);
    vnl_matrix<double> vector = mm.V();
    vnl_vector<double> sigvalue = mm.lambdas();
    
    
    vnl_matrix_fixed<double,3,3> result;
    result.fill(0);
    

        result[0][0]=vector[0][0];
        result[0][1]=vector[0][1];
        result[1][1]=vector[1][1];
        result[1][0]=vector[1][0];
    
    return result;
}

} // end namespace itk
#endif
