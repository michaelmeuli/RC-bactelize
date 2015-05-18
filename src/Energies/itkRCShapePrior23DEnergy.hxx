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

#include "itkRCShapePrior23DEnergy.h"
#include "itkSphereBitmapImageSource.h"
#include "itkRCEnergyBaseClass.h"



namespace itk
{
    template<typename TLabelImage, typename TTemplateImage, typename TEnergyDifference>
    RCShapePrior23DEnergy<TLabelImage, TTemplateImage, TEnergyDifference>
    ::RCShapePrior23DEnergy() {
        m_BackgroundValue = 0;
        m_BackgroundLabel = 0;
        m_OrderOfMoments = 20;
    }
    
    template<typename TLabelImage, typename TTemplateImage, typename TEnergyDifference>
    void RCShapePrior23DEnergy<TLabelImage, TTemplateImage, TEnergyDifference>
    ::PrintSelf(std::ostream & os, Indent indent) const
    {
        Superclass::PrintSelf(os, indent);
        os << indent << "Background value: " << m_BackgroundValue <<  std::endl;
        os << indent << "Background label:" << m_BackgroundLabel << std::endl;
        os << indent << "Order of moments: " << m_OrderOfMoments << std::endl; 
    }
    
    template<typename TLabelImage, typename TTemplateImage, typename TEnergyDifference>
    void RCShapePrior23DEnergy<TLabelImage, TTemplateImage, TEnergyDifference>
    ::PrepareEnergyCalculationForIteration() {
        

        //RenewShapePrior();

        /*
        //new version
        unsigned int vAddSize=m_PointsToBeAdded.size();
        unsigned int vDeleteSize=m_PointsToBeDeleted.size();

        //vClocka=clock();
        vnl_matrix<unsigned int> vAddCoordinate(vAddSize,2);
        vnl_matrix<unsigned int> vDeleteCoordinate(vDeleteSize,2);
        for (unsigned int i=0;i<vAddSize;i++)
        {
         vAddCoordinate(i,0)=m_PointsToBeAdded[i][0];
         vAddCoordinate(i,1)=m_PointsToBeAdded[i][1];
        }
        for (unsigned int i=0;i<vDeleteSize;i++)
        {
         vDeleteCoordinate(i,0)=m_PointsToBeDeleted[i][0];
         vDeleteCoordinate(i,1)=m_PointsToBeDeleted[i][1];
        }
        //vClockb=clock();
        //std::cout<<*m_Moments2D.GetCentralMoments()<<std::endl;
        if (vAddSize != 0 || vDeleteSize !=0)
        {
         m_Moments2D.StepUpdate(&vAddCoordinate, &vDeleteCoordinate);
         m_PointsToBeAdded.clear();
         m_PointsToBeDeleted.clear();
        }
        //std::cout<<*m_Moments2D.GetMoments()<<std::endl;
        //vClockc=clock();
        //std::cout<<"data change time: "<<vClockb-vClocka<<std::endl;
        //std::cout<<"Update Moments: "<<vClockc-vClockb<<std::endl;
         */


        // the old version by scanning
        /*here the code is rescanning the label image*/
        ///*
            
        typename LabelToMoments2DMapType::iterator vActiveLabelsIt =
                m_Moments2Ds.begin();
        typename LabelToMoments2DMapType::iterator vActiveLabelsItEnd =
                m_Moments2Ds.end();         
        for (; vActiveLabelsIt != vActiveLabelsItEnd; ++vActiveLabelsIt) {
            LabelAbsPixelType vLabelOfInterest = vActiveLabelsIt->first;
            int vCount = m_Count[vLabelOfInterest];
            vnl_matrix<unsigned int> vCoordinate;
            vCoordinate.set_size(vCount, 2);
            int vx, vy;
            IndexType vCurrentIndex;
            typedef itk::ImageRegionConstIterator<LabelImageType> LabelImageConstIteraterType;
            LabelImageConstIteraterType vLabelIt(this->m_LabelImage, this->m_LabelImage->GetBufferedRegion());
            vCount = 0;
            clock_t vClocka, vClockb, vClockc; //, vClockd;
            vClocka = clock();
            for (vLabelIt.GoToBegin(); !vLabelIt.IsAtEnd(); ++vLabelIt) {
                LabelPixelType vLabelAtPositionOfTheIterator = abs(vLabelIt.Get());
                if (vLabelAtPositionOfTheIterator == vLabelOfInterest) {
                    vCurrentIndex = vLabelIt.GetIndex();
                    vx = vCurrentIndex[0];
                    vy = vCurrentIndex[1];
                    vCoordinate(vCount, 0) = vx;
                    vCoordinate(vCount, 1) = vy;
                    vCount++;
                }
            }
            
            unsigned int vWidth = this->m_LabelImage->GetLargestPossibleRegion().GetSize()[0];
            unsigned int vHight = this->m_LabelImage->GetLargestPossibleRegion().GetSize()[1];
            //for 2D moments, set the depth 1
            unsigned int vDepth = 1;
            //std::cout<<"**** current points number ***"<<vCount<<std::endl;
            m_Moments2Ds[vLabelOfInterest]->SetInputData(&vCoordinate, vWidth, vHight, vDepth);
            //                vClockb = clock();
            m_Moments2Ds[vLabelOfInterest]->Legendre();
            //                vClockc = clock();
            //            double gong_test = (*m_Moments2D.GetLegendreMoments()-*m_Moments2D.GetReference()).fro_norm();
            //std::cout<<"**** current moments ***"<<std::endl<<*m_Moments2D.GetLegendreMoments()<<std::endl;
            
            //                 FILE *pf= fopen("/Users/ygong/energy.txt","a");
            //                 fprintf(pf,"%f scanning: %d compute: %d\n",gong_test, (int)(vClockb-vClocka),(int)(vClockc-vClockb));
            //                 fclose(pf);
            //                 
            //std::cout<<" *** distance ***"<<gong_test<<std::endl;
            //std::cout<<"**** reference ***"<<std::endl<<*m_Moments2D.GetReference()<<std::endl;
            //
        }
    }
    
    template<typename TLabelImage, typename TTemplateImage, typename TEnergyDifference>
    void RCShapePrior23DEnergy<TLabelImage, TTemplateImage, TEnergyDifference>
    ::PrepareEnergyCalculation() {
        
        ShapePriorReferenceInit();
        //InitializeShapePrior(this->GetInitInput(),&m_Moments2D);
    }
        
    template<typename TLabelImage, typename TTemplateImage, typename TEnergyDifference>
    inline 
    typename RCShapePrior23DEnergy<TLabelImage, TTemplateImage, TEnergyDifference>::
    InternalEnergyReturnType
    RCShapePrior23DEnergy<TLabelImage, TTemplateImage, TEnergyDifference>
    ::EvaluateEnergyDifference(IndexType aInd, 
            LabelPixelType aLabelFrom, 
            LabelPixelType aLabelTo) {
        
             EnergyDifferenceType vShapeEnergy = 0;


        /*  // Access all the contour points via the m_InnerContourContainer:
          OuterContourContainerIteratorType vPointIterator = m_InnerContourContainer.begin();
          OuterContourContainerIteratorType vPointIteratorEnd = m_InnerContourContainer.end();
		
		
                  int vCount=0;
          for(; vPointIterator != vPointIteratorEnd; ++vPointIterator) {
              ContourIndexType vCurrentIndex = vPointIterator->first; // x,y,(z) position
              //int vXPos = vCurrentIndex[0];
              //OuterContourContainerValueType vContourParticle = vPointIterator->second;
              //LabelAbsPixelType vAbsLabel = vContourParticle.m_label;
		
                          vCount++;
          }
                  vPointIterator = m_InnerContourContainer.begin();
                  vnl_matrix<double> vCoordinate;
                  vnl_matrix<double> vCoordinateBefore(vCount,2);
		
                  int vXPos, vYPos;
                  if (aFrom != 0)
                  {
                          vCoordinate.set_size(vCount-1,2);
                          vCount=0;

                          for(; vPointIterator != vPointIteratorEnd; ++vPointIterator)
                          {
                                  ContourIndexType vCurrentIndex = vPointIterator->first; // x,y,(z) position
                                  vXPos = vCurrentIndex[0];
                                  vYPos = vCurrentIndex[1];
                                  vCoordinateBefore(vCount,0)=vXPos;
                                  vCoordinateBefore(vCount,1)=vYPos;
                                  if (vXPos==aIndex[0]&&vYPos==aIndex[1]) continue;
                                  vCoordinate(vCount,0)=vXPos;
                                  vCoordinate(vCount,1)=vYPos;
                                  vCount++;
                          }
                  }
                  else
                  {
                          for(; vPointIterator != vPointIteratorEnd; ++vPointIterator)
                          {
                                  ContourIndexType vCurrentIndex = vPointIterator->first; // x,y,(z) position
                                  vXPos = vCurrentIndex[0];
                                  vYPos = vCurrentIndex[1];
                                  vCoordinateBefore(vCount,0)=vXPos;
                                  vCoordinateBefore(vCount,1)=vYPos;
                                  vCoordinate(vCount,0)=vXPos;
                                  vCoordinate(vCount,1)=vYPos;
                                  vCount++;
                          }
                          vCoordinate(vCount,0)=vXPos;
                          vCoordinate(vCount,1)=vYPos;
                  }
         */
        /* // Access the label image (via iterator):*/

        /*
        int vCount=m_Count[1];

        //vnl_matrix<double> vCoordinateBefore(vCount,2);
        vnl_matrix<double> vCoordinate;
        vCoordinate.set_size(vCount,2);
        int vx, vy;
        LabelImageIndexType vCurrentIndex;
        typedef itk::ImageRegionConstIterator<LabelImageType> LabelImageConstIteraterType;
             LabelImageConstIteraterType vLabelIt(m_LabelImage, m_LabelImage->GetBufferedRegion());

				
             vCount=0;
             for(vLabelIt.GoToBegin(); !vLabelIt.IsAtEnd(); ++vLabelIt)
             {
                     LabelPixelType vLabelAtPositionOfTheIterator = vLabelIt.Get();
                     if(vLabelAtPositionOfTheIterator == 1)
                     {
                             vCurrentIndex = vLabelIt.GetIndex();
                             vx=vCurrentIndex[0];
                             vy=vCurrentIndex[1];
                             vCoordinate(vCount,0)=vx;
                             vCoordinate(vCount,1)=vy;
                             vCount++;
                     }
             }

            unsigned int vWidth=m_LabelImage->GetLargestPossibleRegion().GetSize()[0];
            unsigned int vHight=m_LabelImage->GetLargestPossibleRegion().GetSize()[1];

            m_Moments2D.SetInputData(&vCoordinate,vWidth,vHight);
            m_Moments2D.Legendre();
            vShapeEnergy=m_Moments2D.Distance(m_Moments2D.GetMoments(), m_Moments2D.GetReference(), false);
         */

        //  
        clock_t vClocka, vClockb;
        vClocka = clock();
            
        if(aLabelTo != 0) {
            vShapeEnergy += m_Moments2Ds[aLabelTo]->AddOnePoint(aInd[0], aInd[1]) * m_Count[aLabelTo];
        }
        if(aLabelFrom != 0) {
            vShapeEnergy += m_Moments2Ds[aLabelFrom]->DeleteOnePoint(aInd[0], aInd[1]) * m_Count[aLabelFrom];
        }
        vClockb = clock();

        /*
                clock_t vClocka, vClockb;
                vClocka = clock();
                if (aTo == 0) {
                    std::cout << "[" << aIndex[0] << "," << aIndex[1] << "]" << std::endl;
                    aMoments.DeleteOnePoint(aIndex[0], aIndex[1]);
                    vShapeEnergy = aMoments.GetDeleteMeasure();
                    //just let it grow
                    //vShapeEnergy=100;
                } else {
                    aMoments.AddOnePoint(aIndex[0], aIndex[1]);
                    vShapeEnergy = aMoments.GetAddMeasure();
                }
                vClockb = clock();
                //std::cout<<"dealwith one point needs time: "<<vClockb-vClocka<<std::endl;
                //vcl_cerr<<aIndex[0]<<" "<<aIndex[1]<<" point energy "<<vShapeEnergy<<std::endl;
         */
        return (this->m_Coefficient * vShapeEnergy);
        

    }
    
    template<typename TLabelImage, typename TTemplateImage, typename TEnergyDifference>
    inline void RCShapePrior23DEnergy<TLabelImage, TTemplateImage, TEnergyDifference>
    ::AddPoint(IndexType aInd, LabelAbsPixelType aLabel) {
        
        if(aLabel == m_BackgroundLabel) return;
        
        CountStatisticsIteratorType vIt = m_Count.find(aLabel);
        if(vIt != m_Count.end()) {
            vIt->second++;
            //update the moments map
            //Gong_test
            //             m_PointsToBeAdded.push_back(aInd);
        } else {
            m_Count[aLabel] = 1;
            
            m_Moments2Ds[aLabel] = Moments2DType::New();
            m_Moments2Ds[aLabel]->SetOrder(m_OrderOfMoments);
            m_Moments2Ds[aLabel]->Init();
        }
        
    }
    
    template<typename TLabelImage, typename TTemplateImage, typename TEnergyDifference>
    inline void RCShapePrior23DEnergy<TLabelImage, TTemplateImage, TEnergyDifference>
    ::RemovePoint(IndexType aInd, LabelAbsPixelType aLabel) {
        
        if(aLabel == m_BackgroundLabel) return;
        assert(m_Count.find(aLabel) != m_Count.end());
        if(--m_Count[aLabel] == 0) {
            KillRegion(aLabel);
        } else {
            //update the moments map
            //Gong_test
            //            m_PointsToBeDeleted.push_back(aInd);

        }
        
    }
    
    template<typename TLabelImage, typename TTemplateImage, typename TEnergyDifference>
    void RCShapePrior23DEnergy<TLabelImage, TTemplateImage, TEnergyDifference>
    ::KillRegion(LabelAbsPixelType aLabel) {
        m_Moments2Ds.erase(aLabel);
        m_Count.erase(aLabel);
    }
    

    template<typename TLabelImage, typename TTemplateImage, typename TEnergyDifference>
    void RCShapePrior23DEnergy<TLabelImage, TTemplateImage, TEnergyDifference>
    ::ShapePriorReferenceInit() {
        
        itk::ImageRegionConstIterator<TemplateImageType> it(m_TemplateImage,
                m_TemplateImage->GetRequestedRegion());
        IndexType vCurrentIndex;
        TemplatePixelType vCurrentPixel;
        int count = 0;
        for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
            vCurrentPixel = it.Get();
            if (vCurrentPixel > (unsigned char) 0) { // jc: gong, what is this cast?
                count++;
            }
        }
        vnl_matrix<unsigned int> ref(count, 2);
        count = 0;
        for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
            vCurrentPixel = it.Get();
            if (vCurrentPixel > (unsigned char) 0) { // jc: gong, what is this cast?
                vCurrentIndex = it.GetIndex();
                ref(count, 0) = vCurrentIndex[0];
                ref(count, 1) = vCurrentIndex[1];
                count++;
            }
        }
        unsigned int vWidth = m_TemplateImage->GetLargestPossibleRegion().GetSize()[0];
        unsigned int vHight = m_TemplateImage->GetLargestPossibleRegion().GetSize()[1];
        //set a vDepth as 1 for 2D
        unsigned int vDepth = 1;

        Moments2DPointerType vMoments2DDummy = Moments2DType::New();
        vMoments2DDummy->SetOrder(m_OrderOfMoments);
        vMoments2DDummy->SetInputData(&ref, vWidth, vHight, vDepth);
        vMoments2DDummy->Init();
        vMoments2DDummy->Legendre();
        vMoments2DDummy->SetReference(vMoments2DDummy->GetMoments());
//        vMoments2DDummy->
        
        CountStatisticsIteratorType vActiveLabelsIt =
                m_Count.begin();
        CountStatisticsIteratorType vActiveLabelsItEnd =
                m_Count.end();         
        for (; vActiveLabelsIt != vActiveLabelsItEnd; ++vActiveLabelsIt) {
            LabelAbsPixelType vLabelOfInterest = vActiveLabelsIt->first;
            m_Moments2Ds[vLabelOfInterest]->SetReference(vMoments2DDummy->GetMoments());
        }
        
    }

    
//    template<typename TLabelImage, typename TTemplateImage, typename TEnergyDifference>
//    void RCShapePrior23DEnergy<TLabelImage, TTemplateImage, TEnergyDifference>
//    ::InitializeShapePrior(InitImageType * aInitImage, Moments2D * aMoments) {
//        itk::ImageRegionIterator<InitImageType> it(aInitImage, aInitImage->GetRequestedRegion());
//        InitImageIndexType vCurrentIndex;
//        InitPixelType vCurrentPixel;
//        int count = 0;
//        for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
//            vCurrentPixel = it.Get();
//            if (vCurrentPixel > (unsigned char) 100) {
//                count++;
//            }
//        }
//        vnl_matrix<unsigned int> ref(count, 2);
//        count = 0;
//        for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
//            vCurrentPixel = it.Get();
//            if (vCurrentPixel > (unsigned char) 100) {
//                vCurrentIndex = it.GetIndex();
//                ref(count, 0) = vCurrentIndex[0];
//                ref(count, 1) = vCurrentIndex[1];
//                count++;
//            }
//
//        }
//        unsigned int vWidth = aInitImage->GetLargestPossibleRegion().GetSize()[0];
//        unsigned int vHight = aInitImage->GetLargestPossibleRegion().GetSize()[1];
//        //set a vDepth as 1 for 2D
//        unsigned int vDepth = 1;
//
//        (*aMoments).SetInputData(&ref, vWidth, vHight, vDepth);
//        (*aMoments).Legendre();
//        //std::cout<<"Set Reference "<<*m_Moments2D.GetReference()<<std::endl;hello
//    }
    
    template<typename TLabelImage, typename TTemplateImage, typename TEnergyDifference>
    typename RCShapePrior23DEnergy<TLabelImage, TTemplateImage, TEnergyDifference>
    ::InternalEnergyReturnType 
    RCShapePrior23DEnergy<TLabelImage, TTemplateImage, TEnergyDifference>
    ::CalcFullEnergy(){
                /*
                if (m_UseShapePrior) {
                    int vCount = m_Count[1];
                    //vcl_cerr<<"The number of region points: "<<vCount<<std::endl;
                    //vnl_matrix<double> vCoordinateBefore(vCount,2);
                    vnl_matrix<double> vCoordinate;
                    vCoordinate.set_size(vCount, 2);
                    int vx, vy;
                    LabelImageIndexType vCurrentIndex;
                    typedef itk::ImageRegionConstIterator<LabelImageType> LabelImageConstIteraterType;
                    LabelImageConstIteraterType vLabelIt(m_LabelImage, m_LabelImage->GetBufferedRegion());


                    vCount = 0;
                    for (vLabelIt.GoToBegin(); !vLabelIt.IsAtEnd(); ++vLabelIt) {
                        LabelPixelType vLabelAtPositionOfTheIterator = abs(vLabelIt.Get());
                        if (vLabelAtPositionOfTheIterator == 1) {
                            vCurrentIndex = vLabelIt.GetIndex();
                            vx = vCurrentIndex[0];
                            vy = vCurrentIndex[1];
                            vCoordinate(vCount, 0) = vx;
                            vCoordinate(vCount, 1) = vy;
                            vCount++;
                        }
                    }
                    unsigned int vWidth = m_LabelImage->GetLargestPossibleRegion().GetSize()[0];
                    unsigned int vHight = m_LabelImage->GetLargestPossibleRegion().GetSize()[1];

                    m_Moments2D.SetInputData(&vCoordinate, vWidth, vHight);
                    m_Moments2D.ComputeLegendreMoments();
                    vEnergy += (*m_Moments2D.GetLegendreMoments()-*m_Moments2D.GetReference()).fro_norm();


                }
         */
        return 0;
    }
        
} // end namespace itk