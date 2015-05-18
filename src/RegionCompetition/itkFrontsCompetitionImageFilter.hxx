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
#ifndef __ITKFRONTSCOMPETITIONIMAGEFILTER_CXX_
#define __ITKFRONTSCOMPETITIONIMAGEFILTER_CXX_


#include "itkFrontsCompetitionImageFilter.h"


namespace itk {

    template <class TInputImage, class TInitImage, class TOutputImage >
    FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    FrontsCompetitionImageFilter() {

        /**
         * This filter needs a data and an initialization image as input.
         */
        this->SetNumberOfRequiredInputs(2);

        /**
         * Initialize topology related members
         */
        m_NeighborhoodSize_FG_Connectivity =
                ForegroundConnectivityType::GetInstance().GetNumberOfNeighbors();
        m_NeighborsOffsets_FG_Connectivity =
                ForegroundConnectivityType::GetInstance().GetNeighborsITKOffsets();

        m_NeighborhoodSize_BG_Connectivity =
                BackgroundConnectivityType::GetInstance().GetNumberOfNeighbors();
        m_NeighborsOffsets_BG_Connectivity =
                BackgroundConnectivityType::GetInstance().GetNeighborsITKOffsets();

        m_TopologicalNumberFunction =
                TopologicalNumberImageFunction<LabelImageType, ForegroundConnectivityType>::New();
        m_TopologicalNumberFunction->SetComputeBackgroundTN(false);


        /**
         * Initialize control members
         */
        m_MaxNbIterations = 100;
        m_converged = false;
        m_spatiallyConverged = true;

        m_AllowFusion = true;
        m_AllowFission = true;
        m_AllowHandles = true;

        m_UseForbiddenRegion = false;
        m_UseFastEvolution = false;
        m_UseSobolevGradients = false;


        m_ForbiddenRegionLabel = NumericTraits<LabelPixelType>::max();

        for (unsigned int vD = 0; vD < m_Dim; vD++) {
            m_SobolevKernelSigma[vD] = 10;
        }
        m_iteration_counter = 0;

        m_AcceptedPointsFactor = 1;
        m_AcceptedPointsReductionFactor = 0.5;

        m_MemoryLength = 10;
       
        
        m_GPUEnergy = NULL;
        m_MinimalCellSize.Fill(10);
    }

    template <class TInputImage, class TInitImage, class TOutputImage >
    FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    ~FrontsCompetitionImageFilter() {
    }

    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    PrintSelf(std::ostream& aOS, Indent aIndent) const {

        Superclass::PrintSelf(aOS, aIndent);
        aOS << aIndent << "Maximal number of iterations: " << m_MaxNbIterations << std::endl;
        aOS << aIndent << "Allow fusion: " << m_AllowFusion << std::endl;
        aOS << aIndent << "Allow fission: " << m_AllowFission << std::endl;
        aOS << aIndent << "Allow handles: " << m_AllowHandles << std::endl;
        aOS << aIndent << "Memory length: " << m_MemoryLength << std::endl;
        aOS << aIndent << "Sobolev kernel sigma: " << m_SobolevKernelSigma << std::endl;
        aOS << aIndent << "Use Sobolev gradients: " << m_UseSobolevGradients << std::endl;
    }

    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    GenerateData() {
        /*
         * Set up the regions and allocate the output-image
         */
        CopyRegionsAndAllocateOutput();

        /*
         * Set up the memory image
         */
        m_MemoryImage = MemoryImageType::New();
        m_MemoryImage->SetRegions(this->GetDataInput()->GetLargestPossibleRegion());
        m_MemoryImage->Allocate();
        m_MemoryImage->FillBuffer(0);
        m_MemoryImage->SetSpacing(this->GetDataInput()->GetSpacing());

        /*
         * Set up the containers and allocate the label image:
         */
        InitializeLabelImageAndContourContainer(this->GetInitInput());
        
        /*
         * Set the inputs for the image and energy functions:
         */
        m_TopologicalNumberFunction->SetInputImage(m_LabelImage);

        for(unsigned int vI = 0; vI < m_ExtEnergies.size(); vI++) {
            m_ExtEnergies[vI]->SetLabelImage(m_LabelImage);
            m_ExtEnergies[vI]->SetDataImage(this->GetInput());
        }
        for(unsigned int vI = 0; vI < m_IntEnergies.size(); vI++) {
            m_IntEnergies[vI]->SetLabelImage(m_LabelImage);
        }
        for(unsigned int vI = 0; vI < m_PIEnergies.size(); vI++) {
            m_PIEnergies[vI]->SetLabelImage(m_LabelImage);
        }
        if(m_GPUEnergy.IsNotNull()) {
            m_GPUEnergy->SetLabelImage(m_LabelImage);
            m_GPUEnergy->SetDataImage(this->GetInput());
        }
        
        /*
         * Initialize standard statistics (mean, variances, length, area etc)
         */
        RenewStatistics(m_LabelImage);
        

        /*
         * Depending on the functional to use, prepare stuff for faster computation.
         */
        PrepareEnergyCaluclation();

//        	WriteInputImageTypeToFile("DEBUG_data.tif", this->GetDataInput(), 1000);
        //      WriteInputImageTypeToFile("DEBUG_psf.tif",m_PSF,100000); 

        m_iteration_counter = 0;

        


        bool vConvergence = false;
        int vCountdownToAnnealing = m_MemoryLength / 2 + 1;
        if (this->GetDebug()) {
            WriteDebugImages(m_iteration_counter);
        }

        /*
         * Main loop of the RegionCompetition
         */
        while (m_MaxNbIterations > m_iteration_counter && !(vConvergence)) {
            m_iteration_counter++;
            std::cout << ".";
            std::flush(std::cout);

            m_spatiallyConverged = true;

            itk::TimeProbe vIterationTimer;
            vIterationTimer.Start();

            vConvergence = DoOneIteration();

            vIterationTimer.Stop();

            m_TimingsHist.push_back(vIterationTimer.GetTotal());

            if (m_spatiallyConverged == true) {
                vCountdownToAnnealing -= 1; // enter the annealing phase
                std::cout << "Spatially converged." << std::endl;
            }
            if (vCountdownToAnnealing < 0) {
                m_AcceptedPointsFactor *= m_AcceptedPointsReductionFactor;
                std::cout << "Annealing: nb of accepted points reduced to: "
                        << m_AcceptedPointsFactor << std::endl;
                vCountdownToAnnealing = m_MemoryLength / 2 + 1;
            }

            if (this->GetDebug()) {
                std::cout << "iteration: " << m_iteration_counter << std::endl;
                std::cout << "---------------" << std::endl;
                WriteDebugImages(m_iteration_counter);
            }
        }
        m_converged = vConvergence;



        /*
         * Write the label image in a convenient form to the filters output image
         */
        WriteLabelImageToOutput();

//        std::ofstream vTimingsOutStream;
//        vTimingsOutStream.open("history_timings.txt");
//        for (unsigned int i = 0; i < m_TimingsHist.size(); i++) {
//            vTimingsOutStream << i + 1 << "\t" << m_TimingsHist[i] << "\n";
//        }
//        vTimingsOutStream.close();



        for(unsigned int vI = 0; vI < m_ExtEnergies.size(); vI++) {
            m_ExtEnergies[vI]->CleanUp();
        }
        for(unsigned int vI = 0; vI < m_IntEnergies.size(); vI++) {
            m_IntEnergies[vI]->CleanUp();   
        }
        for(unsigned int vI = 0; vI < m_PIEnergies.size(); vI++) {
            m_PIEnergies[vI]->CleanUp();  
        }
        if(m_GPUEnergy.IsNotNull()) {
            m_GPUEnergy->CleanUp();
        }

        if (m_converged) {
            std::cout << "convergence after " << m_iteration_counter <<
                    "iterations." << std::endl;
        } else {
            std::cout << "no convergence !" << std::endl;
        }
    }


    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    WriteContourParticleContainer(const char* aFilename, CellListedHashMapType &aContainer) {
        std::ofstream vOutStream;
        vOutStream.open(aFilename);
        CellListedHashMapIteratorType vIt = aContainer.begin();
        CellListedHashMapIteratorType vItEnd = aContainer.end();
        for (; vIt != vItEnd; ++vIt) {
            vOutStream << vIt->first << "\t" << vIt->second.m_label << "\t"
                    << vIt->second.m_candidateLabel << "\t" << vIt->second.m_energyDifference
                    << "\t" << vIt->second.m_modifiedCounter << "\n";
        }
        vOutStream.close();
    }

    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    DebugCellListedHashMapIterator(const char* aFilename, CellListedHashMapType &aContainer) {
        std::ofstream vOutStream;
        vOutStream.open(aFilename);

        //        typename OuterContourContainerType::const_iterator vIt = aContainer.begin();
        //        typename OuterContourContainerType::const_iterator vItEnd = aContainer.end();
        CellListedHashMapIteratorType vIt = aContainer.begin();
        CellListedHashMapIteratorType vItEnd = aContainer.end();
        for (; vIt != vItEnd; ++vIt) {
            CellListedHashMapCellIteratorType vNeighIt = m_InnerContourContainer.cell_begin(vIt->first);
            CellListedHashMapCellIteratorType vNeighEnd = m_InnerContourContainer.cell_end(vIt->first);
            for (; vNeighIt != vNeighEnd; ++vNeighIt) {
                vOutStream << vIt->first << "\t" << vIt.getCellID() << "\t" << vNeighIt->first << "\t" << vNeighIt.getCellID() << "\n";
            }
        }
        vOutStream.close();
    }

    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    WriteDebugImages(unsigned int aIterationNumber) {
//        std::stringstream vSS;
//        vSS << "debugLabelsAtIt_" << (aIterationNumber + 10000) << ".tif";
//        std::string vS = vSS.str();
//        WriteLabelImageToFile(const_cast<char *> (vS.c_str()));
        
        std::stringstream vSSDebug2;
        vSSDebug2 << "debugLabelsAtIt2_" << (aIterationNumber + 10000) << ".tif";
        WriteLabelImageToFile2(const_cast<char *> (vSSDebug2.str().c_str()));

//        std::stringstream vSSMem;
//        vSSMem << "memoryAtIt_" << (aIterationNumber + 10000) << ".tif";
//        std::string vSMem = vSSMem.str();
//        WriteMemoryImageToFile(const_cast<char *> (vSMem.c_str()));

    }

    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    WriteEnergyDifferencesToImageFile(
    CellListedHashMapType * aContainer, const char* aFilename) {
        typedef unsigned short RGBComponentType;
        typedef RGBPixel<RGBComponentType> PxType;
        typedef itk::Image<PxType, m_Dim> ImType;
        typename ImType::Pointer vGradientImage = ImType::New();

        vGradientImage->SetRegions(m_LabelImage->GetLargestPossibleRegion());
        vGradientImage->Allocate();
        vGradientImage->SetSpacing(this->GetDataInput()->GetSpacing());
        PxType vBlackRGBPixel;
        vBlackRGBPixel.Set(0, 0, 0);
        vGradientImage->FillBuffer(vBlackRGBPixel);

        //iterate through the sample to be able to scale to the output image range
        typename ContourContainerValueType::EnergyDifferenceType vEMax = 0;
        CellListedHashMapIteratorType vCIt = aContainer->begin();
        for (; vCIt != aContainer->end(); ++vCIt) {
            if (fabs(vCIt->second.m_energyDifference) > vEMax) {
                vEMax = fabs(vCIt->second.m_energyDifference);
            }
        }

        vCIt = aContainer->begin();
        for (; vCIt != aContainer->end(); ++vCIt) {

            RGBComponentType vScaledEnergyDiffValue = static_cast<RGBComponentType> (
                    fabs(vCIt->second.m_energyDifference) /
                    vEMax * itk::NumericTraits<RGBComponentType>::max());

            //            RGBComponentType vScaledEnergyDiffValue = static_cast<RGBComponentType> (
            //                    fabs(vCIt->second.m_energyDifference));

            //            RGBComponentType vScaledEnergyDiffValue = fabs(vCIt->second.m_energyDifference);

            PxType vRGBPixel;
            vRGBPixel.Set(0, 0, 0);
            if (vCIt->second.m_energyDifference > 0) {
                vRGBPixel.SetGreen(vScaledEnergyDiffValue);
            } else {
                vRGBPixel.SetRed(vScaledEnergyDiffValue);
            }
            vGradientImage->SetPixel(vCIt->first, vRGBPixel);
        }

        std::stringstream vSS;
        vSS << "energyDifferencesAtIt_" << (m_iteration_counter + 10000) << ".tif";

        typedef ImageFileWriter<ImType> FileWriterType;
        typename FileWriterType::Pointer vFileWriter = FileWriterType::New();
        vFileWriter->SetInput(vGradientImage);
        vFileWriter->SetFileName(aFilename);
        vFileWriter->Update();
    }

    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    WriteMaskImageTypeToFile(const char* aFilename, MaskImageType* aImage, float aScale) {
        typedef short LocalOutputPixelType;
        typedef itk::Image<LocalOutputPixelType, m_Dim> LocalOutputImageType;

        typedef ShiftScaleImageFilter<MaskImageType, MaskImageType> ShiftScaleFilterType;
        typename ShiftScaleFilterType::Pointer vShiftScaleFilter = ShiftScaleFilterType::New();
        vShiftScaleFilter->SetInput(aImage);
        vShiftScaleFilter->SetScale(aScale);

        typedef CastImageFilter<MaskImageType, LocalOutputImageType> CasterType;
        typename CasterType::Pointer vCaster = CasterType::New();
        vCaster->SetInput(vShiftScaleFilter->GetOutput());

        typedef ImageFileWriter<LocalOutputImageType> FileWriterType;
        typename FileWriterType::Pointer vFileWriter = FileWriterType::New();
        vFileWriter->SetInput(vCaster->GetOutput());
        vFileWriter->SetFileName(aFilename);
        vFileWriter->Update();
    }

    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    WriteInputImageTypeToFile(const char* aFilename, InputImageType* aImage, float aScale) {
        typedef short LocalOutputPixelType;
        typedef itk::Image<LocalOutputPixelType, m_Dim> LocalOutputImageType;

        typedef ShiftScaleImageFilter<InternalImageType, InternalImageType> ShiftScaleFilterType;
        typename ShiftScaleFilterType::Pointer vShiftScaleFilter = ShiftScaleFilterType::New();
        vShiftScaleFilter->SetInput(aImage);
        vShiftScaleFilter->SetScale(aScale);

        typedef CastImageFilter<InternalImageType, LocalOutputImageType> CasterType;
        typename CasterType::Pointer vCaster = CasterType::New();
        vCaster->SetInput(vShiftScaleFilter->GetOutput());

        typedef ImageFileWriter<LocalOutputImageType> FileWriterType;
        typename FileWriterType::Pointer vFileWriter = FileWriterType::New();
        vFileWriter->SetInput(vCaster->GetOutput());
        vFileWriter->SetFileName(aFilename);
        vFileWriter->Update();
    }

    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    WriteLabelImageTypeToFile(const char* aFilename, LabelImageType* aImage) {

        typedef short LocalOutputPixelType;
        typedef itk::Image<LocalOutputPixelType, m_Dim> LocalOutputImageType;

        typedef CastImageFilter<LabelImageType, LocalOutputImageType> CasterType;
        typename CasterType::Pointer vCaster = CasterType::New();
        vCaster->SetInput(aImage);

        typedef ImageFileWriter<LocalOutputImageType> FileWriterType;
        typename FileWriterType::Pointer vFileWriter = FileWriterType::New();
        vFileWriter->SetInput(vCaster->GetOutput());
        vFileWriter->SetFileName(aFilename);
        vFileWriter->Update();
    }

    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    WriteLabelAbsImageTypeToFile(const char* aFilename, LabelAbsImageType* aImage) {

        typedef unsigned short LocalOutputPixelType;
        typedef itk::Image<LocalOutputPixelType, m_Dim> LocalOutputImageType;

        typedef CastImageFilter<LabelAbsImageType, LocalOutputImageType> CasterType;
        typename CasterType::Pointer vCaster = CasterType::New();
        vCaster->SetInput(aImage);

        typedef ImageFileWriter<LocalOutputImageType> FileWriterType;
        typename FileWriterType::Pointer vFileWriter = FileWriterType::New();
        vFileWriter->SetInput(vCaster->GetOutput());
        vFileWriter->SetFileName(aFilename);
        vFileWriter->Update();
    }

    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    WriteMemoryImageToFile(char* aFilename) {
        typedef unsigned char OutputPixelType2;
        typedef itk::Image<OutputPixelType2, m_Dim> OutputImageType2;

        typedef CastImageFilter<MemoryImageType, OutputImageType2> CasterType;
        typename CasterType::Pointer vCaster = CasterType::New();

        typedef ImageFileWriter<OutputImageType2> FileWriterType;
        typename FileWriterType::Pointer vFileWriterSink = FileWriterType::New();

        vCaster->SetInput(m_MemoryImage);

        vFileWriterSink->SetInput(vCaster->GetOutput());
        vFileWriterSink->SetFileName(aFilename);
        vFileWriterSink->SetUseCompression(false);
        vFileWriterSink->Update();
    }

    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    WriteLabelImageToFile(char* aFilename) {

        typedef unsigned short OutputPixelType2;
        typedef itk::Image<OutputPixelType2, m_Dim> OutputImageType2;

        typedef CastImageFilter<LabelImageType, OutputImageType2> CasterType;
        typename CasterType::Pointer vCaster = CasterType::New();

        typedef ImageFileWriter<OutputImageType2> FileWriterType;
        typename FileWriterType::Pointer vFileWriterSink = FileWriterType::New();

        vCaster->SetInput(m_LabelImage);
        vCaster->Update();

        CellListedHashMapIteratorType vCIt = m_InnerContourContainer.begin();
        for (; vCIt != m_InnerContourContainer.end(); ++vCIt) {
            vCaster->GetOutput()->SetPixel(vCIt->first, vCIt->second.m_label);
            //            vCaster->GetOutput()->SetPixel(vCIt->first, vCIt->second.m_label + 500);
        }

        vFileWriterSink->SetInput(vCaster->GetOutput());
        vFileWriterSink->SetFileName(aFilename);
        vFileWriterSink->SetUseCompression(false);
        vFileWriterSink->Update();
    }

    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    WriteLabelImageToFile2(char* aFilename) {

        typedef unsigned short DebugPixelType;
        typedef itk::Image<DebugPixelType, m_Dim> DebugImageType;


        typedef ImageFileWriter<DebugImageType> FileWriterType;
        typename FileWriterType::Pointer vFileWriterSink = FileWriterType::New();

        typename DebugImageType::Pointer vOutImage = DebugImageType::New();
        vOutImage->SetRegions(m_LabelImage->GetBufferedRegion());
        vOutImage->Allocate();
        vOutImage->SetSpacing(this->GetDataInput()->GetSpacing());

        ImageRegionIterator<DebugImageType> vOutputIt(vOutImage,
                vOutImage->GetBufferedRegion());
        ImageRegionConstIterator<LabelImageType> vLabelIt(m_LabelImage,
                m_LabelImage->GetBufferedRegion());


        for (vOutputIt.GoToBegin(), vLabelIt.GoToBegin(); !vOutputIt.IsAtEnd();
                ++vLabelIt, ++vOutputIt) {
            if (vLabelIt.Get() < 0) {
                //                vOutputIt.Set(static_cast<OutputPixelType>(m_NLabels+10));
                vOutputIt.Set(static_cast<OutputPixelType> (2));
            } else if (vLabelIt.Get() == 0) {
                vOutputIt.Set(static_cast<OutputPixelType> (0));
            } else if (vLabelIt.Get() == m_ForbiddenRegionLabel) {
                vOutputIt.Set(static_cast<OutputPixelType> (3));
            } else {
                vOutputIt.Set(static_cast<OutputPixelType> (1));
            }
        }
        vFileWriterSink->SetInput(vOutImage);
        vFileWriterSink->SetFileName(aFilename);
        vFileWriterSink->SetUseCompression(false);
        vFileWriterSink->Update();
    }
    
 

    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    UpdateMemoryImage() {
        /// let the memory forget...
        ImageRegionIterator<MemoryImageType> vMemIt(m_MemoryImage,
                m_MemoryImage->GetBufferedRegion());
        for (vMemIt.GoToBegin(); !vMemIt.IsAtEnd(); ++vMemIt) {
            MemoryPixelType vV = vMemIt.Get();
            if (vV > 0) {
                vMemIt.Set(vV - 1);
            }
        }
    }
    
     
    template <class TInputImage, class TInitImage, class TOutputImage >
    bool FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    DoOneIteration() {
        PrepareEnergyCaluclationForEachIteration();
//        RemoveSinglePointRegions();
        UpdateMemoryImage();
        bool vConvergenceA = EvaluateCandidatesAndApplyMoves();
//        CleanUp();
        return vConvergenceA;
    }

    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    WriteLabelImageToOutput() {
        ImageRegionIterator<TOutputImage> vOutputIt(this->GetOutput(),
                this->GetOutput()->GetBufferedRegion());
        ImageRegionConstIterator<LabelImageType> vLabelIt(m_LabelImage,
                this->GetInput()->GetBufferedRegion());

        for (vOutputIt.GoToBegin(), vLabelIt.GoToBegin(); !vOutputIt.IsAtEnd();
                ++vLabelIt, ++vOutputIt) {
            if (vLabelIt.Get() == m_ForbiddenRegionLabel) {
                vOutputIt.Set(static_cast<OutputPixelType> (0));
                continue;
            }
            vOutputIt.Set(static_cast<OutputPixelType> (abs(vLabelIt.Get())));
            //            vOutputIt.Set(static_cast<OutputPixelType> (
            //                    255*m_Means[abs(vLabelIt.Get())]));
        }
    }




    /**
     * AddNeighborsAtRemove adds all the necessary neighbors to the contour container
     * when removing a point from the region and adding it to the background
     * region.
     * Basically all the FG-Neighborhood neighbors with the same label are added
     * to the container.
     */
    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    AddNeighborsAtRemove(
    LabelAbsPixelType aAbsLabel,
    ContourContainerKeyType aIndex) {

        LabelImageNeighborhoodIteratorRadiusType vLabelImageIteratorRadius;
        vLabelImageIteratorRadius.Fill(1);

        LabelImageNeighborhoodIteratorType vLabelImageIterator(vLabelImageIteratorRadius,
                m_LabelImage,
                m_LabelImage->GetBufferedRegion());

        vLabelImageIterator.SetLocation(aIndex);

        for (unsigned int vN = 0; vN < m_NeighborhoodSize_FG_Connectivity; vN++) {
            InputImageOffsetType vOff = m_NeighborsOffsets_FG_Connectivity[vN];
            LabelPixelType vL = vLabelImageIterator.GetPixel(vOff);
            if (vL > 0 && static_cast<unsigned int> (vL) == aAbsLabel) {
                ContourIndexType vI = aIndex;
                vI += vOff;
                ContourContainerValueType vAddV;
                vAddV.m_label = aAbsLabel;
                vAddV.m_candidateLabel = 0;
                //                vAddV.m_modifiedCounter = 0;
                vAddV.m_data = this->GetInput()->GetPixel(vI);
                //                (*aTempAddConntainer)[vI] = vAddV;

                m_LabelImage->SetPixel(vI, -aAbsLabel);
                m_InnerContourContainer[vI] = vAddV;
            } else if (vL < 0 && static_cast<unsigned int> (vL) == aAbsLabel) {
                /// the point is already in the contour. We reactivate it by
                /// ensuring that the energy is calculated in the next iteration:
                ContourIndexType vI = aIndex;
                vI += vOff;
                m_InnerContourContainer[vI].m_modifiedCounter = 0;
            }
        }
    }

    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    UpdateStatisticsWhenJump(
    InternalImageIndexType aIndex,
    InputPixelType aData,
    LabelAbsPixelType aFromLabel,
    LabelAbsPixelType aToLabel) {

        InternalPixelType vCurrentImageValue = aData;
        
        for(unsigned int vI = 0; vI < m_ExtEnergies.size(); vI++) {
            m_ExtEnergies[vI]->SwitchPoint(aIndex, aFromLabel, aToLabel, vCurrentImageValue);
        }
        for(unsigned int vI = 0; vI < m_IntEnergies.size(); vI++) {
            m_IntEnergies[vI]->SwitchPoint(aIndex, aFromLabel, aToLabel);             
        }
        for(unsigned int vI = 0; vI < m_PIEnergies.size(); vI++) {
            m_PIEnergies[vI]->SwitchPoint(aIndex, aFromLabel, aToLabel);             
        }
        if(m_GPUEnergy.IsNotNull()) {
            m_GPUEnergy->SwitchPoint(aIndex, aFromLabel, aToLabel, vCurrentImageValue);    
        }

    }

 
    /**
     * Maintain the inner contour container:
     * - Remove all the indices in the BG-connected neighborhood, that are
     *   interior points, from the contour container.
     *   Interior here means that none neighbors in the FG-Neighborhood
     *   has a different label.
     */
    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    MaintainNeighborsAtAdd(
    LabelAbsPixelType aLabelAbs,
    ContourContainerKeyType aIndex) {

        LabelImageNeighborhoodIteratorRadiusType vLabelImageIteratorRadius;
        vLabelImageIteratorRadius.Fill(1);

        // we set the pixel value already to ensure the that the 'enclosed' check
        // afterwards works.
        m_LabelImage->SetPixel(aIndex, -aLabelAbs);

        for (unsigned int vI = 0; vI < m_NeighborhoodSize_BG_Connectivity; vI++) {

            InputImageOffsetType vOff = m_NeighborsOffsets_BG_Connectivity[vI];
            aIndex += vOff;

            // we use the label image here to identify the old points in the
            // inner container. This early-abort hopefully speeds-up the execution (instead of
            // just deleting all FG-neighbors from the container).
            //            if (vLabelImageIterator.GetCenterPixel() == -static_cast<int> (aLabelAbs) &&
            //                    IsEnclosedByLabel(vLabelImageIterator)) {
            if (m_LabelImage->GetPixel(aIndex) == -static_cast<int> (aLabelAbs) &&
                    IsEnclosedByLabel_FGConnectivity(aIndex, aLabelAbs)) {
                /// It might happen that a point, that was already accepted as a 
                /// candidate, gets enclosed by other candidates. This
                /// points must not be added to the container afterwards and thus
                /// removed from the main list.
                m_InnerContourContainer.erase(aIndex);
                //                aCandidatatesContourContainer->erase(aIndex);
                //                vLabelImageIterator.SetCenterPixel(aLabelAbs);
                m_LabelImage->SetPixel(aIndex, aLabelAbs);
                //                m_LabelImage->SetPixel(vCurrentIndex, vLabelAbs); //TODO: use the iterator above (make it non-const)
            }

            aIndex -= vOff;
        }

        if (IsEnclosedByLabel_FGConnectivity(aIndex, aLabelAbs)) {
            m_InnerContourContainer.erase(aIndex);
            m_LabelImage->SetPixel(aIndex, aLabelAbs);
        }

    }


    template <class TInputImage, class TInitImage, class TOutputImage >
    typename FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    EnergyDifferenceType
    FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    CalculateAverageInternalMotherEnergyForLabel(
    CellListedHashMapType* aCandidateContainer,
    ContourContainerValueType* aParticle,
    unsigned int aLabel) {
        EnergyDifferenceType vAvgInnerEnergies = static_cast<EnergyDifferenceType> (0);
        EnergyDifferenceType vNbMothersOfPropLabel = static_cast<EnergyDifferenceType> (0);

        typedef typename ContourParticleType::MotherIndexListType::iterator
        MotherListIteratorType;
        MotherListIteratorType vMothersIt = aParticle->m_motherIndices.begin();
        MotherListIteratorType vMothersItEnd = aParticle->m_motherIndices.begin();

        for (; vMothersIt != vMothersItEnd; ++vMothersIt) {
            ContourContainerValueType& vMotherRef = aCandidateContainer->find(*vMothersIt)->second;
            if (vMotherRef.m_label == aLabel) {
                vAvgInnerEnergies += vMotherRef.m_energyDifference;
                vNbMothersOfPropLabel++;
            }
        }

        return vAvgInnerEnergies / vNbMothersOfPropLabel;
    }
    

    template <class TInputImage, class TInitImage, class TOutputImage >
    typename FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    EnergyDifferenceType
    FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    CalculateInternalEnergyDifferenceForLabel(
    ContourContainerKeyType aContourIndex,
    ContourContainerValueType * aContourParticlePtr,
    LabelAbsPixelType aToLabel) {

        EnergyDifferenceType vEnergy = 0;
        LabelAbsPixelType vCurrentLabel = aContourParticlePtr->m_label;

        for(unsigned int vI = 0; vI < m_PIEnergies.size(); vI++) {
            // TODO: the following distinction should happen in the energies and
            // not here. 
            if(aToLabel != 0) { // daughter
                vEnergy += CalculateAverageInternalMotherEnergyForLabel(
                        &m_InnerContourContainer, aContourParticlePtr, aToLabel);  
            } else {
                vEnergy += m_PIEnergies[vI]->EvaluateEnergyDifference(aContourIndex, 
                        m_InnerContourContainer.cell_begin(aContourIndex),
                        m_InnerContourContainer.cell_end(aContourIndex));
            }
        }
        for(unsigned int vI = 0; vI < m_IntEnergies.size(); vI++) {
            vEnergy += m_IntEnergies[vI]->EvaluateEnergyDifference(aContourIndex, vCurrentLabel, aToLabel);            
        }

        return vEnergy;
    }



    template <class TInputImage, class TInitImage, class TOutputImage >
    typename FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    EnergyDifferenceAndMergePairType
    FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    CalculateExternalEnergyDifferenceForLabel(
    ContourContainerKeyType aContourIndex,
    ContourContainerValueType * aContourParticlePtr,
    LabelAbsPixelType aToLabel) {

        EnergyDifferenceType vEnergy = 0;
        EnergyEventType vMerge = false;

        InputPixelType vCurrentImageValue = aContourParticlePtr->m_data;
        LabelAbsPixelType vCurrentLabel = aContourParticlePtr->m_label;

        /// Calculate the change in energy due to the change of intensity when changing
        /// from one label 'from' to another 'to'.
        for(unsigned int vI = 0; vI < m_ExtEnergies.size(); vI++) {
            typename ExternalEnergyType::ExternalEnergyReturnType vV = m_ExtEnergies[vI]->
            EvaluateEnergyDifference(aContourIndex,vCurrentLabel, aToLabel, vCurrentImageValue);
            
            vEnergy += vV.first;
            vMerge = vV.second || vMerge;
        }

        return EnergyDifferenceAndMergePairType(vEnergy, vMerge);
    }

    template <class TInputImage, class TInitImage, class TOutputImage >
    float FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    SumAllEnergies(CellListedHashMapType* aContainer) {
        CellListedHashMapIteratorType vPointIterator = aContainer->begin();
        CellListedHashMapIteratorType vPointIteratorEnd = aContainer->end();
        EnergyDifferenceType vTotalEnergyDiff = 0;
        for (vPointIterator = aContainer->begin(); vPointIterator != vPointIteratorEnd; ++vPointIterator) {
            vTotalEnergyDiff += vPointIterator->second.m_energyDifference;
        }
        return vTotalEnergyDiff;
    }

    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    FilterCandidatesContainerUsingRanks(CellListedHashMapType* aContainer) {

        if (aContainer->size() > 0) {
            // Copy the candidates to a set (of ContourParticleWithIndex). This
            // will sort them according to their energy gradients.
            typedef ContourParticleWithIndex<ContourContainerKeyType, ContourContainerValueType>
                    ContourParticleWithIndexType;
            typedef std::list<ContourParticleWithIndexType> ContourParticleWithIndexListType;
            ContourParticleWithIndexListType vSortedList;
            CellListedHashMapIteratorType vPointIterator = aContainer->begin();
            CellListedHashMapIteratorType vPointIteratorEnd = aContainer->end();

            for (vPointIterator = aContainer->begin(); vPointIterator != vPointIteratorEnd; ++vPointIterator) {
                ContourParticleWithIndexType vCand(vPointIterator->first, vPointIterator->second);
                vSortedList.push_back(vCand);
            }

            vSortedList.sort();

            unsigned int vNbElements = vSortedList.size();
            vNbElements = static_cast<int> (static_cast<float> (vNbElements)
                    * m_AcceptedPointsFactor + 0.5);

            /// Fill the container with the best candidate first, then
            /// the next best that does not intersect the tabu region of
            /// all inserted points before.
            aContainer->clear();
            typename ContourParticleWithIndexListType::iterator vSortedListIterator = vSortedList.begin();
            for (; vNbElements >= 1 && vSortedListIterator != vSortedList.end(); ++vSortedListIterator) {
                vNbElements--;
                ContourContainerKeyType vCandCIndex = vSortedListIterator->m_ContourIndex;
                CellListedHashMapIteratorType vAcceptedCandIterator = aContainer->begin();
                bool vValid = true;
                for (; vAcceptedCandIterator != aContainer->end(); ++vAcceptedCandIterator) {
                    ContourContainerKeyType vCIndex = vAcceptedCandIterator->first;
                    
                }
                if (vValid) {
                    /// This candidate passed the test and is added to the TempRemoveCotainer:
                    (*aContainer)[vSortedListIterator->m_ContourIndex] = vSortedListIterator->m_ContourParticle;
                }
            }
        }
    }


    /**
     * The method performs all steps needed to switch a pixel to a different
     * label. No tests about topology or validity of the arguments are performed.
     * Concretely,
     * - the label image is updated.
     * - the main contour container is maintained.
     * - Statistics are updated, including energy specific tasks (maintaining the
     *   deconvolution-model image.
     * aNewCCPointsContainer contains contour points that have been added to the
     *
     */
    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    ChangeContourParticleLabelToCandidateLabel(ContourIndexType aIndex, ContourContainerValueType* aParticle) {

        LabelAbsPixelType vFromLabel = aParticle->m_label;
        LabelAbsPixelType vToLabel = aParticle->m_candidateLabel;

        ///
        /// The particle was modified,reset the counter in order to process
        /// the particle in the next iteration.
        ///
        aParticle->m_modifiedCounter = 0;

        ///
        /// Update the label image. The new point is either a contour point or 0,
        /// therefor the negative label value is set.
        ///
        m_LabelImage->SetPixel(aIndex, -static_cast<LabelPixelType> (vToLabel));

        ///
        /// Update the statistics of the propagating and the loser region.
        ///
        UpdateStatisticsWhenJump(aIndex, aParticle->m_data, vFromLabel, vToLabel);

        ///
        /// SPATIAL CONVERGENCE MANAGEMENT
        ///
        if (m_MemoryImage->GetPixel(aIndex) == 0) {
            /// We did not recently visit this index, so the front is still moving
            m_spatiallyConverged = false;
        }
        m_MemoryImage->SetPixel(aIndex, m_MemoryLength);


        /// TODO: A bit a dirty hack: we store the old label for the relabeling
        ///       procedure later on...either introduce a new variable or rename the
        ///       variable (which doesn't work currently :-).
        aParticle->m_candidateLabel = vFromLabel;

        ///
        /// Clean up
        ///

        /// The loser region (if it is not the BG region) has to add the
        /// neighbors of the lost point to the contour list.
        if (vFromLabel != 0) {
            AddNeighborsAtRemove(vFromLabel, aIndex);
        }

        ///
        /// Erase the point from the surface container in case it now belongs to
        /// the background. Else, add the point to the container (or replace it
        /// in case it has been there already).
        ///
        if (vToLabel == 0) {
            m_InnerContourContainer.erase(aIndex);
        } else {
            ContourContainerValueType vContourParticle = *aParticle;
            vContourParticle.m_label = vToLabel;
            /// The point may or may not exist already in the m_InnerContainer.
            /// The old value, if it exists, is just overwritten with the new
            /// contour point (with a new label).
            m_InnerContourContainer[aIndex] = vContourParticle;
        }

        /// Remove 'enclosed' contour points from the container. For the BG this
        /// makes no sense.
        if (vToLabel != 0) {
            MaintainNeighborsAtAdd(vToLabel, aIndex);
        }
    }

    //    /**
    //     * RemoveContourParticle performs all steps to remove a contour point (switch
    //     * to the BG label). It adds  the new ContourParticles to aNewContourParticles container.
    //     * (points that are now adjacent to another label).
    //     */
    //    template <class TInputImage, class TInitImage, class TOutputImage >
    //    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    //    RemoveContourParticle(
    //    OuterContourContainerType* aNewContourParticles,
    //    LabelAbsPixelType aFromLabel,
    //    OuterContourContainerIteratorType aIt) {
    //        OuterContourContainerKeyType vCurrentIndex = aIt->first;
    //        UpdateStatisticsWhenJump(aIt, aFromLabel, 0);
    //        m_LabelImage->SetPixel(vCurrentIndex, 0);
    //        AddNeighborsAtRemove(aNewContourParticles, aFromLabel, vCurrentIndex);
    //        m_InnerContourContainer.erase(vCurrentIndex);
    //
    //        if (m_EnergyFunctional == e_DeconvolutionPC) {
    //            UpdateConvolvedImage(vCurrentIndex, aFromLabel, 0);
    //        }
    //    }
    //
    //    template <class TInputImage, class TInitImage, class TOutputImage >
    //    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    //    AddContourParticle(
    //    OuterContourContainerType* aTempAddContainer,
    //    OuterContourContainerIteratorType aIt) {
    //
    //        OuterContourContainerKeyType vCurrentIndex = aIt->first;
    //        (*aTempAddContainer)[vCurrentIndex] = aIt->second; //add it
    //        LabelAbsPixelType vToLabel = aIt->second.m_candidateLabel;
    //
    //        /// The update of the label image here is necessary to detect non-simple
    //        /// points of 2 'approaching fronts' (fronts that are exactly 2 pixel
    //        /// apart from each other and both regions 'grow').
    //        m_LabelImage->SetPixel(vCurrentIndex, -vToLabel);
    //        UpdateStatisticsWhenJump(aIt, 0, vToLabel);
    //        MaintainNeighborsAtAdd(vToLabel, vCurrentIndex, aTempAddContainer);
    //
    //        if (m_EnergyFunctional == e_DeconvolutionPC) {
    //            UpdateConvolvedImage(vCurrentIndex, 0, vToLabel);
    //        }
    //    }

    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    AddUpEnergies(CellListedHashMapType* aCandidateContainer) {

        /// smooth the external energy terms 
        if (m_UseSobolevGradients) {

            /// This will write the filtered values in the filtereEnergyDifference
            /// field of the particles:
            SmoothExternalEnergies(aCandidateContainer);

            /// Add up the external and internal energy terms
            CellListedHashMapIteratorType vPointIterator = aCandidateContainer->begin();
            CellListedHashMapIteratorType vPointsEnd = aCandidateContainer->end();
            for (; vPointIterator != vPointsEnd; ++vPointIterator) {
                vPointIterator->second.m_energyDifference +=
                        vPointIterator->second.m_filteredEnergyDifference;
            }
        } else {
            CellListedHashMapIteratorType vPointIterator = aCandidateContainer->begin();
            CellListedHashMapIteratorType vPointsEnd = aCandidateContainer->end();
            for (; vPointIterator != vPointsEnd; ++vPointIterator) {
                vPointIterator->second.m_energyDifference +=
                        vPointIterator->second.m_externalEnergyDifference;
            }
        }
    }

    template <class TInputImage, class TInitImage, class TOutputImage >
    bool FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    EvaluateCandidatesAndApplyMoves() {
        /// Convergence is set to false if a point moved:
        bool vConvergence = true;

        /// Clear the old candidate set
        m_Candidates.clear();
        m_Seeds.clear();


        /// clear the competing regions map, it will be refilled in
        /// RebuildCandidateList:
        m_CompetingRegionsMap.clear();

        /// Build the new candidate set, this includes computation of energies:
        RebuildAndOptimizeCandidateList(&m_Candidates);

        /// Sum up the (maybe smoothed) external and internal energies.
        AddUpEnergies(&m_Candidates);
        

        /// Debugging output
        unsigned int vNumberAllCandidates = m_Candidates.size();
        unsigned int vNumberMovedCandidates = 0;

        if (this->GetDebug()) {
//            std::stringstream vSS2;
//            vSS2 << "ContourParticlesAfterRebuildAtIt_" << (m_iteration_counter + 10000) << ".txt";
//            WriteContourParticleContainer(vSS2.str().c_str(), m_Candidates);

            //            std::stringstream vSS4;
            //            vSS4 << "cellListedHashMap_it_" << (m_iteration_counter + 10000) << ".txt";
            //            DebugCellListedHashMapIterator(vSS4.str().c_str(), m_AllCandidates);

//            std::stringstream vSS3;
//            vSS3 << "ContourParticles_" << (m_iteration_counter + 10000) << ".txt";
            //            WriteContourParticleContainer(vSS3.str().c_str(), m_InnerContourContainer);

//            std::stringstream vSS;
//            vSS << "energyDifferencesAtIt_" << (m_iteration_counter + 10000) << ".tif";
//            WriteEnergyDifferencesToImageFile(&m_Candidates, vSS.str().c_str());
//
//            std::cout << "nb candidates: " << m_Candidates.size() << std::endl;
        }

        /**
         * Find topologically compatible candidates and store their indices in
         * vLegalIndices.
         */
        itk::TimeProbe vTimerFrontPropagation;
        vTimerFrontPropagation.Start();


        typedef std::list<ContourIndexType> ContourIndexListType;
        ContourIndexListType vLegalIndices;
        ContourIndexListType vIllegalIndices;
        typedef ContourParticleWithIndex<ContourContainerKeyType, ContourContainerValueType>
                ContourParticleWithIndexType;
        typedef std::list<ContourParticleWithIndexType> NetworkMembersListType;

        CellListedHashMapIteratorType vPointIterator = m_Candidates.begin();
        CellListedHashMapIteratorType vPointsEnd = m_Candidates.end();
        for (; vPointIterator != vPointsEnd; ++vPointIterator) {

            /// Check if this point was processed already
            if (!vPointIterator->second.m_processed) {

                /// Check if it is a mother: only mothers can be seed points
                /// of topological networks. Daughters are always part of a
                /// topo network of a mother.
                if (!vPointIterator->second.m_isMother)
                    continue;

                /**
                 * Build the dependency network for this seed point:
                 */
                std::stack<ContourIndexType> vIndicesToVisit;
                NetworkMembersListType vSortedNetworkMembers;
                vIndicesToVisit.push(vPointIterator->first);
                vPointIterator->second.m_processed = true;

                while (!vIndicesToVisit.empty()) {

                    ContourIndexType vSeedIndex = vIndicesToVisit.top();
                    vIndicesToVisit.pop();

//                    std::cout << "building network for " << vSeedIndex << std::endl;
                    ContourContainerValueType vCurrentMother = m_Candidates.find(vSeedIndex)->second;

                    /// Add the seed point to the network
                    ContourParticleWithIndexType vSeedContourParticleWithIndex(
                            vSeedIndex, vCurrentMother);
                    vSortedNetworkMembers.push_back(vSeedContourParticleWithIndex);

                    // Iterate all children of the seed, push to the stack if there
                    // is a mother.
                    typename ContourParticleType::DaughterIndexListType::iterator vDaughterIt =
                            vCurrentMother.m_daughterIndices.begin();
                    typename ContourParticleType::DaughterIndexListType::iterator vDaughterItEnd =
                            vCurrentMother.m_daughterIndices.end();

                    for (; vDaughterIt != vDaughterItEnd; ++vDaughterIt) {
                        ContourIndexType vDaughterContourIndex = *vDaughterIt;

//                        ContourContainerValueType& vDaughterContourParticle =
//                                m_Candidates[vDaughterContourIndex];
//                        std::cout << "daughter is " << vDaughterContourIndex << std::endl;
                        ContourContainerValueType& vDaughterContourParticle =
                                m_Candidates.find(vDaughterContourIndex)->second;
//                        std::cout << " with associated particle: " << vDaughterContourParticle << std::endl; 
                        if (!vDaughterContourParticle.m_processed) {
                            vDaughterContourParticle.m_processed = true;
                            
                            if (vDaughterContourParticle.m_isMother) {
                                vIndicesToVisit.push(vDaughterContourIndex);
                            } else {
                                ContourParticleWithIndexType vDaughterContourParticleWithIndex(
                                        vDaughterContourIndex, vDaughterContourParticle);
                                vSortedNetworkMembers.push_back(vDaughterContourParticleWithIndex);
                            }

                            /// Push all the non-processed mothers of this daughter to the stack
                            typename ContourParticleType::MotherIndexListType::iterator vDMIt =
                                    vDaughterContourParticle.m_motherIndices.begin();
                            typename ContourParticleType::MotherIndexListType::iterator vDMItEnd =
                                    vDaughterContourParticle.m_motherIndices.end();
                            for (; vDMIt != vDMItEnd; ++vDMIt) {
                                ContourIndexType vDM = *vDMIt;

//                                ContourContainerValueType& vMotherOfDaughterPoint = m_Candidates[vDM];
                                ContourContainerValueType& vMotherOfDaughterPoint = m_Candidates.find(vDM)->second;
                                if (!vMotherOfDaughterPoint.m_processed) {
                                    vMotherOfDaughterPoint.m_processed = true;
                                    vIndicesToVisit.push(vDM);
                                }
                            }
                        }
                    }
                }

                /**
                 * sort the network
                 */
                vSortedNetworkMembers.sort();

                /**
                 * Filtering: Accept all members in ascending order that are
                 * compatible with the already selected members in the network.
                 */
                typedef std::set< ContourIndexType > ContourIndexSetType;
                ContourIndexSetType vSelectedCandidateIndices;

                typename NetworkMembersListType::iterator vNetworkIt = vSortedNetworkMembers.begin();
                typename NetworkMembersListType::iterator vNetworkItEnd = vSortedNetworkMembers.end();


                //                ContourIndexType vFirstIndex = vNetworkIt->m_ContourIndex;
                //                vSelectedCandidateIndices.insert(vFirstIndex);
                //                //                ContourIndexType vDummyIndex;
                //                //                vSelectedCandidateIndices.insert(vDummyIndex);
                //                ++vNetworkIt;
                for (; vNetworkIt != vNetworkItEnd; ++vNetworkIt) {
                    /// If a mother is accepted, the reference count of all the
                    /// daughters (with the same label) has to be decreased.
                    /// Rules: a candidate in the network is a legal candidate if:
                    /// - If (daughter): The reference count >= 1. (Except the
                    ///                  the candidate label is the BG - this allows
                    ///                  creating BG regions inbetween two competing
                    ///                  regions).
                    /// - If ( mother ): All daughters (with the same 'old' label) in the
                    ///   accepted list have still a reference count > 1.
                    bool vLegalMove = true;
                    
                    
                    ///
                    /// RULE 1: If c is a daughter point, the reference count r_c is > 0.
                    ///
                    if (vNetworkIt->m_ContourParticle.m_isDaughter) {
                        ContourContainerValueType& vCand = 
                                m_Candidates.find(vNetworkIt->m_ContourIndex)->second;
                        if (vCand.m_referenceCount < 1 && vCand.m_candidateLabel != 0) {
                            vLegalMove = false;
                        }
                    }

                    ///
                    /// RULE 2: All daughters already accepted the label of this
                    ///        mother have at least one another mother.
                    /// AND
                    /// RULE 3: Mothers are still valid mothers (to not introduce
                    ///         holes in the FG region).
                    ///
                    if (vLegalMove && vNetworkIt->m_ContourParticle.m_isMother) {
                        /// Iterate the daughters and check their reference count
                        typename ContourParticleWithIndexType::ContourParticleType::MotherIndexListType::iterator
                        vDaughterIndicesIterator = vNetworkIt->m_ContourParticle.m_daughterIndices.begin();
                        typename ContourParticleWithIndexType::ContourParticleType::MotherIndexListType::iterator
                        vDaughterIndicesIteratorEnd = vNetworkIt->m_ContourParticle.m_daughterIndices.end();

                        bool vRule3Fulfilled = false;

                        for (; vDaughterIndicesIterator != vDaughterIndicesIteratorEnd;
                                ++vDaughterIndicesIterator) {

//                            ContourContainerValueType& vDaughterPoint =
//                                    m_Candidates[*vDaughterIndicesIterator];
                            ContourContainerValueType& vDaughterPoint =
                                    m_Candidates.find(*vDaughterIndicesIterator)->second;

                            /// rule 2:
                            typename ContourIndexSetType::iterator vAcceptedDaugtherIt =
                                    vSelectedCandidateIndices.find(*vDaughterIndicesIterator);

                            if (vAcceptedDaugtherIt != vSelectedCandidateIndices.end()) {
                                /// This daughter has been accepted and needs
                                /// a reference count > 1, else the move is
                                /// invalid.
                                if (vDaughterPoint.m_candidateLabel == vNetworkIt->m_ContourParticle.m_label &&
                                        vDaughterPoint.m_referenceCount <= 1) {
                                    vLegalMove = false;
                                    break;
                                }
                            }

                            ///
                            /// rule 3:
                            if (!vRule3Fulfilled) {

                                if (vAcceptedDaugtherIt == vSelectedCandidateIndices.end()) {
                                    /// There is a daughter that has not yet been accepted.
                                    vRule3Fulfilled = true;
                                } else {
                                    /// the daughter has been accepted, but may
                                    /// have another candidate label(rule 3b):                                       
                                    if (m_Candidates.find(*vAcceptedDaugtherIt)->second.m_candidateLabel !=
                                            vNetworkIt->m_ContourParticle.m_label) {
                                        vRule3Fulfilled = true;
                                    }
                                }
                            }
                        }

                        if (!vRule3Fulfilled) vLegalMove = false;
                    }

                    if (vLegalMove) {
                        /// the move is legal, store the index in the container
                        /// with accepted candidates of this network.
                        vSelectedCandidateIndices.insert(vNetworkIt->m_ContourIndex);

                        /// also store the candidate in the global candidate
                        /// index container:
                        vLegalIndices.push_back(vNetworkIt->m_ContourIndex);

                        /// decrease the references of its daughters(with the same 'old' label).
                        typename ContourParticleWithIndexType::ContourParticleType::MotherIndexListType::iterator
                        vDaughterIndicesIterator = vNetworkIt->m_ContourParticle.m_daughterIndices.begin();
                        typename ContourParticleWithIndexType::ContourParticleType::MotherIndexListType::iterator
                        vDaughterIndicesIteratorEnd = vNetworkIt->m_ContourParticle.m_daughterIndices.end();
                        for (; vDaughterIndicesIterator != vDaughterIndicesIteratorEnd;
                                ++vDaughterIndicesIterator) {
                            ContourContainerValueType& vDaughterPoint =
                                    m_Candidates.find(*vDaughterIndicesIterator)->second;
                            if (vDaughterPoint.m_candidateLabel ==
                                    vNetworkIt->m_ContourParticle.m_label) {
                                vDaughterPoint.m_referenceCount--;
                            }
                        }
                    } else {
                        vIllegalIndices.push_back(vNetworkIt->m_ContourIndex);
                    }
                }
            }
        }
        vTimerFrontPropagation.Stop();
        m_TimeHistFrontPropagation.push_back(vTimerFrontPropagation.GetTotal());


        /**
         * Filter all candidates with the illigal indices
         */
        typename ContourIndexListType::iterator vIlligalIndicesIt = vIllegalIndices.begin();
        typename ContourIndexListType::iterator vIlligalIndicesItEnd = vIllegalIndices.end();
        for (; vIlligalIndicesIt != vIlligalIndicesItEnd; ++vIlligalIndicesIt) {
            m_Candidates.erase(*vIlligalIndicesIt);
        }

        /**
         * Filter candidates according to their energy
         */
        ContourIndexListType vIndicesWithPosEnergyDiff;
        vPointIterator = m_Candidates.begin();
        vPointsEnd = m_Candidates.end();
        for (;vPointIterator != vPointsEnd; ++vPointIterator) {
//            CellListedHashMapIteratorType vStoreIt = vPointIterator; // iterator to work with
//            ++vPointIterator; // safely increment (before erasing anything in the container)
//            if (vStoreIt->second.m_energyDifference >= 0) {
//                m_Candidates.erase(vStoreIt);
//            }
            if(vPointIterator->second.m_energyDifference >= 0) {
                vIndicesWithPosEnergyDiff.push_back(vPointIterator->first);
            }
        }
        typename ContourIndexListType::iterator vPosEgyIt = vIndicesWithPosEnergyDiff.begin();
        for(;vPosEgyIt != vIndicesWithPosEnergyDiff.end(); ++vPosEgyIt){
            m_Candidates.erase(*vPosEgyIt);
        }
        

        /**
         * Intermediate step: filter the candidates according to their rank
         * and spacial position.
         */
        if (m_AcceptedPointsFactor < 0.99) {
            FilterCandidatesContainerUsingRanks(&m_Candidates);
        }

        /**
         * Move all the points that are simple. Non simple points remain in the
         * candidates list.
         */
        typedef typename TopologicalNumberCalculatorType::ForegroundTopologicalNumbersType
        FGTNType;
        typedef typename FGTNType::iterator FGTNIteratorType;
        FGTNType vFGTNvector;
        bool vChange = true;
        /// We first move all the FG-simple points. This we do because it happens
        /// that points that are not simple at the first place get simple after
        /// the change of other points. The non-simple points will be treated
        /// in a separate loop afterwards.
        while (vChange && !m_Candidates.empty()) {
            vChange = false;
            vPointIterator = m_Candidates.begin();
            vPointsEnd = m_Candidates.end();

//            while (vPointIterator != vPointsEnd) {
            ContourIndexListType vIndicesToDel;
            for (;vPointIterator != vPointsEnd; ++vPointIterator) {

                ContourContainerKeyType vCurrentIndex = vPointIterator->first;

                vFGTNvector = (m_TopologicalNumberFunction->
                        EvaluateAdjacentRegionsFGTNAtIndex(vCurrentIndex));


                FGTNIteratorType vTopoNbItr;
                FGTNIteratorType vTopoNbItrEnd = vFGTNvector.end();
                bool vSimple = true;
                /// Check for FG-simplicity:
                for (vTopoNbItr = vFGTNvector.begin();
                        vTopoNbItr != vTopoNbItrEnd; ++vTopoNbItr) {
                    if (vTopoNbItr->second.first != 1 || vTopoNbItr->second.second != 1) {
                        // This is a FG simple point; perform the move.
                        vSimple = false;
                    }
                }
                if (vSimple) {
                    vChange = true;
                    ChangeContourParticleLabelToCandidateLabel(vPointIterator->first, &(vPointIterator->second));
                    vNumberMovedCandidates++;
                    vIndicesToDel.push_back(vPointIterator->first);
                    vConvergence = false;
                }

                /// we will reuse the processed flag to indicate if a particle is a seed
                vPointIterator->second.m_processed = false;
            }
            typename ContourIndexListType::iterator vIndToDelIt = vIndicesToDel.begin();
            for(;vIndToDelIt != vIndicesToDel.end(); ++vIndToDelIt){
                m_Candidates.erase(*vIndToDelIt);
            }
        }

        /// Now we know that all the points in the list are 'currently' not simple.
        /// We move them anyway (if topologicial constraints allow) but record
        /// (for every particle) where to relabel (using the seed set). Placing
        /// the seed is necessary for every particle to ensure relabeling even
        /// if a bunch of neighboring particles change. The seed will be ignored
        /// later on if the corresponding FG region is not present in the
        /// neighborhood anymore. 

        /// TODO: The follwoing code is dependent on the iteration order if splits/handles
        /// are not allowed. A solution would be to sort the candidates beforehands.
        /// This should be computationally not too expensive since we assume there
        /// are not many non-simple points.
        //        for (unsigned int vRep = 0; vRep < 2; vRep++) {
        vPointIterator = m_Candidates.begin();
        vPointsEnd = m_Candidates.end();
        for (;vPointIterator != vPointsEnd; ++vPointIterator) {

//            CellListedHashMapIteratorType vStoreIt = vPointIterator; // iterator to work with
//            ++vPointIterator; // safely increment (before erasing anything in the container)

            ContourContainerKeyType vCurrentIndex = vPointIterator->first;
            LabelAbsPixelType vCurrentLabel = vPointIterator->second.m_label;
            LabelAbsPixelType vCandidateLabel = vPointIterator->second.m_candidateLabel;

            bool vValidPoint = true;

            vFGTNvector = (m_TopologicalNumberFunction->
                    EvaluateAdjacentRegionsFGTNAtIndex(vCurrentIndex));

            FGTNIteratorType vTopoNbItr;
            FGTNIteratorType vTopoNbItrEnd = vFGTNvector.end();

            /// Check for handles:
            /// if the point was not disqualified already and we disallow
            /// introducing handles (not only self fusion!), we check if
            /// there is an introduction of a handle.
            if (vValidPoint && !m_AllowHandles) {
                for (vTopoNbItr = vFGTNvector.begin();
                        vTopoNbItr != vTopoNbItrEnd; ++vTopoNbItr) {

                    if (vTopoNbItr->first == vCandidateLabel) {
                        if (vTopoNbItr->second.first > 1) {
                            vValidPoint = false;
                        }
                    }

                    /// criterion to detect surface points or surface junctions.
                    /// or surface-curve junctions. Changing a surface-curve
                    /// junction will also cause a split.

                    if (vTopoNbItr->second.second > 1) {
                        vValidPoint = false;
                    }
                }
            }

            /// Check for splits:
            /// This we have to do either to forbid
            /// the change in topology or to register the seed point for
            /// relabelling.
            /// if the point was not disqualified already and we disallow
            /// splits, then we check if the 'old' label undergoes a split.
            if (vValidPoint) {
                // - "allow introducing holes": T_FG(x, L = l') > 1
                // - "allow splits": T_FG > 2 && T_BG == 1
                bool vSplit = false;
                for (vTopoNbItr = vFGTNvector.begin();
                        vTopoNbItr != vTopoNbItrEnd; ++vTopoNbItr) {
                    if (vTopoNbItr->first == vCurrentLabel) {
                        if (vTopoNbItr->second.first > 1) {
                            vSplit = true;
                        }
                    }
                }
                if (vSplit) {
                    if (m_AllowFission) {
                        RegisterSeedsAfterSplit(m_LabelImage, vCurrentIndex, vCurrentLabel, &m_Candidates);
                    } else {
                        /// disallow the move.
                        vValidPoint = false;
                    }
                }
            }

            /// If the move doesn't change topology or is allowed (and registered
            /// as seed) to change the topology, perform the move (in the
            /// second iteration; in the first iteration seed points need to
            /// be collected):
            if (vValidPoint) {
                ChangeContourParticleLabelToCandidateLabel(vPointIterator->first, &(vPointIterator->second));
                vNumberMovedCandidates++;
                vConvergence = false;

                /// check if the point was a seed and if so, hand the scapegoat
                /// to the another particle
                if (vPointIterator->second.m_processed) {
                    RegisterSeedsAfterSplit(m_LabelImage, vCurrentIndex, vCurrentLabel, &m_Candidates);
                    m_Seeds.erase(SeedType(vCurrentIndex, vCurrentLabel));
                }
            }
        }


        bool vSplit = false;
        bool vMerge = false;
        /// Perform relabeling of the regions that did a split:
        SeedSetIteratorType vSeedIt = m_Seeds.begin();
        SeedSetIteratorType vSeedItEnd = m_Seeds.end();

        unsigned int vNSplits = 0;
        itk::TimeProbe vTimerRelabeling;
        vTimerRelabeling.Start();
        for (; vSeedIt != vSeedItEnd; ++vSeedIt) {
            RelabelRegionsAfterSplit(m_LabelImage, vSeedIt->first, vSeedIt->second);

            vSplit = true;
            vNSplits++;
        }

        /// Merge competing regions if they meet the merging criterion.
        if (m_AllowFusion) {
            typename CompetingRegionsMapType::iterator vCRit = m_CompetingRegionsMap.begin();
            typename CompetingRegionsMapType::iterator vCRend = m_CompetingRegionsMap.end();
            for (; vCRit != vCRend; ++vCRit) {
                LabelAbsPixelType vLabel1 = (vCRit->first)[0];
                LabelAbsPixelType vLabel2 = (vCRit->first)[1];
                m_MergingHist.push_back(m_iteration_counter);
                m_MergingHist.push_back(vLabel1);
                m_MergingHist.push_back(vLabel2);

                RelabelRegionsAfterFusions(m_LabelImage,
                        vCRit->second, vLabel1);
                vMerge = true;
            }

        }
        vTimerRelabeling.Stop();
        m_TimeHistRelabeling.push_back(vTimerRelabeling.GetTotal());

        std::cout << "Moved particles: " << vNumberMovedCandidates << " out of " <<
                vNumberAllCandidates << ": " <<
                100. * vNumberMovedCandidates / vNumberAllCandidates << "%. Timings:" << std::endl;
        std::cout << "Reb." << std::setprecision(2) << m_TimeHistRebuilding.back() << "s\t";
        std::cout << "Egy." << std::setprecision(2) << m_TimeHistEnergyComputation.back() << "\t";

        std::cout << "Min." << std::setprecision(2) << m_TimeHistMinimization.back() << "s\t";
        std::cout << "FrP." << std::setprecision(2) << m_TimeHistFrontPropagation.back() << "s\t";
        std::cout << "Rel." << std::setprecision(2) << m_TimeHistRelabeling.back() << "s\t";
        std::cout << std::endl;

        return vConvergence;
    }

    template <class TInputImage, class TInitImage, class TOutputImage >
    inline bool FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    IsBoundaryPoint(const LabelImageIndexType& aIndex) const {
        LabelAbsPixelType vLabelAbs = static_cast<unsigned int> (abs(m_LabelImage->GetPixel(aIndex)));
        for (unsigned int vI = 0; vI < m_NeighborhoodSize_FG_Connectivity; vI++) {
            if (static_cast<unsigned int> (
                    abs(m_LabelImage->GetPixel(aIndex + m_NeighborsOffsets_FG_Connectivity[vI])))
                    != vLabelAbs) {
                return true;
            }
        }
        return false;
    }

    template <class TInputImage, class TInitImage, class TOutputImage >
    inline bool FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    HasBGNeighbor(const LabelImageIndexType& aIndex) const {

        for (unsigned int vI = 0; vI < m_NeighborhoodSize_FG_Connectivity; vI++) {
            if (m_LabelImage->GetPixel(aIndex + m_NeighborsOffsets_FG_Connectivity[vI]) == 0) {
                return true;
            }
        }
        return false;
    }


        
    template <class TInputImage, class TInitImage, class TOutputImage >
    inline bool FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    IsEnclosedByLabel_FGConnectivity(const LabelImageNeighborhoodIteratorType& aIt) const {

        LabelPixelType vAbsLabel = abs(aIt.GetCenterPixel());

        for (unsigned int vI = 0; vI < m_NeighborhoodSize_FG_Connectivity; vI++) {
            if (abs(aIt.GetPixel(m_NeighborsOffsets_FG_Connectivity[vI])) != vAbsLabel)
                return false;
        }
        return true;
    }

    template <class TInputImage, class TInitImage, class TOutputImage >
    inline bool FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    IsEnclosedByLabel_FGConnectivity(const LabelImageIndexType& aIndex, LabelAbsPixelType aAbsLabel) const {
        for (unsigned int vI = 0; vI < m_NeighborhoodSize_FG_Connectivity; vI++) {
            if (static_cast<unsigned int> (abs(m_LabelImage->GetPixel(
                    aIndex + m_NeighborsOffsets_FG_Connectivity[vI]))) != aAbsLabel)
                return false;
        }
        return true;
    }

    
    


    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    SmoothExternalEnergies(CellListedHashMapType* aContainer) {

        CellListedHashMapIteratorType vPointIterator = aContainer->begin();
        CellListedHashMapIteratorType vPointsEnd = aContainer->end();

        /// caluclate the normalization constant for the multivariate Gaussian
        float vDeterminant = 1;
        for (unsigned int vD = 0; vD < m_Dim; vD++) {
            vDeterminant *= m_SobolevKernelSigma[vD];
        }
        float vIsotropeNormalPDFNormalizer = pow(2.0f * M_PI, -static_cast<float> (m_Dim) / 2.0f) *
                (1.0f / pow(vDeterminant, 0.5));


        for (; vPointIterator != vPointsEnd; ++vPointIterator) {

            CellListedHashMapCellIteratorType vNeighIt = aContainer->cell_begin(vPointIterator->first);
            CellListedHashMapCellIteratorType vNeighItEnd = aContainer->cell_end(vPointIterator->first);

            EnergyDifferenceType vInsideSum = 0;
            EnergyDifferenceType vOutsideSum = 0;

            unsigned int vNbInsideNeighbors = 0;
            unsigned int vNbOutsideNeighbors = 0;
            for (; vNeighIt != vNeighItEnd; ++vNeighIt) {
                if ((vNeighIt->second.m_label == vPointIterator->second.m_label &&
                        vNeighIt->second.m_candidateLabel == vPointIterator->second.m_candidateLabel)
                        ||
                        (vNeighIt->second.m_candidateLabel == vPointIterator->second.m_label &&
                        vNeighIt->second.m_label == vPointIterator->second.m_candidateLabel)) {


                    // calculate the exponent of the isotrope multivariate normal.
                    float vExponent = 0;
                    for (unsigned int vD = 0; vD < m_Dim; vD++) {
                        vExponent += ((vPointIterator->first)[vD] - (vNeighIt->first)[vD]) *
                                ((vPointIterator->first)[vD] - (vNeighIt->first)[vD]) /
                                m_SobolevKernelSigma[vD];
                    }

                    // Put a cutoff (radius)
                    if (vExponent > 9.0f) { // corresponds to a rel error of 0.01.
                        continue;
                    }

                    if (vNeighIt->second.m_label != vPointIterator->second.m_label) {
                        vNbOutsideNeighbors++;
                        vOutsideSum += vNeighIt->second.m_externalEnergyDifference *
                                vIsotropeNormalPDFNormalizer * exp(-vExponent / 2.0f);
                    } else {
                        vNbInsideNeighbors++;
                        vInsideSum += vNeighIt->second.m_externalEnergyDifference *
                                vIsotropeNormalPDFNormalizer * exp(-vExponent / 2.0f);
                    }
                }
            }

            if (vNbInsideNeighbors == 0) vNbInsideNeighbors = 1;
            if (vNbOutsideNeighbors == 0) vNbOutsideNeighbors = 1;

            vPointIterator->second.m_filteredEnergyDifference =
                    vInsideSum / static_cast<float> (vNbInsideNeighbors) -
                    vOutsideSum / static_cast<float> (vNbOutsideNeighbors);

        }
    }



    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    DistroyGPUArrays(
    unsigned int **aGPUCoordinates,
    unsigned int* aGPUToLabel,
    float* aGPUOutInternalEnergy,
    float* aGPUOutExternalEnergy,
    unsigned int* aGPUOutMerge) {
        /// clean up GPU arrays
        for (int i = 0; i < m_Dim; i++)
            delete [] aGPUCoordinates[i];

        delete [] aGPUCoordinates;
        delete [] aGPUToLabel;
        delete [] aGPUOutExternalEnergy;
        delete [] aGPUOutInternalEnergy;
        delete [] aGPUOutMerge;
    }


    template <class TInputImage, class TInitImage, class TOutputImage >
    unsigned int FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    InitializeGPUArrays(
    CellListedHashMapType* aCandidateContainer,
    std::vector<std::vector<GPU_IntType> >* aCoords,
    std::vector<LabelAbsPixelType>* aCandLabels,
    unsigned int ***aGPUCoordinates,
    unsigned int** aGPUToLabel,
    float** aGPUOutInternalEnergy,
    float** aGPUOutExternalEnergy,
    unsigned int** aGPUOutMerge) {

        /// find number of daughters and mothers
        unsigned int vNbOfParticles = aCandLabels->size() + m_InnerContourContainer.size();

        /// pad the number of output arguments to the next multiple of 512:
        vNbOfParticles = vNbOfParticles + (512 - (vNbOfParticles & 511));

        /// Initialize
        //        aGPUCoordinates = new unsigned int*[m_Dim];
        (*aGPUCoordinates) = (unsigned int**) calloc(m_Dim, sizeof (unsigned int*));
        for (unsigned int i = 0; i < m_Dim; i++) {
            unsigned int* aPtr = (unsigned int*) calloc(vNbOfParticles, sizeof (unsigned int));

            (*aGPUCoordinates)[i] = aPtr;
        }

        *aGPUToLabel = new unsigned int[vNbOfParticles];
        *aGPUOutInternalEnergy = new float[vNbOfParticles];
        *aGPUOutExternalEnergy = new float[vNbOfParticles];
        *aGPUOutMerge = new unsigned int[vNbOfParticles];

        for (int i = 0; i < vNbOfParticles; i++) {
            for (unsigned int j = 0; j < m_Dim; j++)
                (*aGPUCoordinates)[j][i] = 0;
            (*aGPUToLabel)[i] = 0;
            (*aGPUOutInternalEnergy)[i] = 0;
            (*aGPUOutExternalEnergy)[i] = 0;
            (*aGPUOutMerge)[i] = 0;
        }

        unsigned int vCandidateCount = 0;

        /// Copy coords and candidate labels of current MOTHERS in cand list
        CellListedHashMapIteratorType vPointIterator = aCandidateContainer->begin();
        CellListedHashMapIteratorType vPointsEnd = aCandidateContainer->end();

        for (; vPointIterator != vPointsEnd; ++vPointIterator) {
            if (vPointIterator->second.m_isMother) {
                ContourContainerKeyType vCurrentIndex = vPointIterator->first;
                for (unsigned int vJ = 0; vJ < m_Dim; vJ++) {
                    (*aGPUCoordinates)[vJ][vCandidateCount] = vCurrentIndex[vJ];
                }
                (*aGPUToLabel)[vCandidateCount] = 0;
                vCandidateCount++;
            }
        }

        /// concat the daughters coords and candidate labels:
        for (unsigned int vI = 0; vI < aCandLabels->size(); vI++) {

            for (unsigned int vJ = 0; vJ < m_Dim; vJ++) {
                (*aGPUCoordinates)[vJ][vCandidateCount] = (*aCoords)[vJ][vI];
            }
            (*aGPUToLabel)[vCandidateCount] = (*aCandLabels)[vI];
            vCandidateCount++;
        }

        return vNbOfParticles;
    }

    
    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    CopyGPUResults(
    CellListedHashMapType* aCandidateContainer,
    std::vector<EnergyDifferenceType>* aInternalEnergy,
    std::vector<EnergyDifferenceType>* aExternalEnergy,
    std::vector<EnergyEventType>* aMerge,
    float* aGPUOutInternalEnergy,
    float* aGPUOutExternalEnergy,
    unsigned int* aGPUOutMerge) {

        /// Copy the values for the mothers:
        CellListedHashMapIteratorType vPointIterator = aCandidateContainer->begin();
        CellListedHashMapIteratorType vPointsEnd = aCandidateContainer->end();
        unsigned int vCandidateCount = 0;
        for (; vPointIterator != vPointsEnd; ++vPointIterator) {
            ContourContainerValueType& vParticle = vPointIterator->second;

            if (vParticle.m_isMother) {
                vParticle.m_energyDifference = aGPUOutInternalEnergy[vCandidateCount];
                vParticle.m_externalEnergyDifference = aGPUOutExternalEnergy[vCandidateCount];

                vCandidateCount++;
            }
        }

        /// Copy the values for the daughters:
        for (unsigned int vI = 0; vI < aInternalEnergy->size(); vCandidateCount++, vI++) {
            (*aMerge)[vI] = aGPUOutMerge[vCandidateCount];
            (*aInternalEnergy)[vI] = aGPUOutInternalEnergy[vCandidateCount];
            (*aExternalEnergy)[vI] = aGPUOutExternalEnergy[vCandidateCount];
        }
    }


    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    AddDaughtersToCandidateList(CellListedHashMapType* aCandidateContainer) {

        /// Iterate the contour list and visit all the neighbors in the
        /// FG-Neighborhood.
        LabelImageNeighborhoodIteratorRadiusType vLabelImageIteratorRadius;
        vLabelImageIteratorRadius.Fill(1);
        LabelImageNeighborhoodIteratorType vLabelImageIterator(vLabelImageIteratorRadius,
                m_LabelImage,
                m_LabelImage->GetBufferedRegion());

        CellListedHashMapIteratorType vPointIterator = m_InnerContourContainer.begin();
        CellListedHashMapIteratorType vPointsEnd = m_InnerContourContainer.end();

        for (; vPointIterator != vPointsEnd; ++vPointIterator) {
            ContourContainerKeyType vCurrentIndex = vPointIterator->first;
            ContourContainerValueType vVal = vPointIterator->second;
            if (m_UseFastEvolution && vVal.m_modifiedCounter - 1 != 0 &&
                    (m_iteration_counter) % (vVal.m_modifiedCounter - 1) != 0) {
                continue;
            }

            vLabelImageIterator.SetLocation(vCurrentIndex);
            LabelAbsPixelType vLabelOfPropagatingRegion = vVal.m_label;
                    
            
            for (unsigned int vI = 0; vI < m_NeighborhoodSize_FG_Connectivity; vI++) {
                InputImageOffsetType vOff = m_NeighborsOffsets_FG_Connectivity[vI];
                LabelAbsPixelType vLabelOfDefender = abs(vLabelImageIterator.GetPixel(vOff));
                if (vLabelOfDefender == static_cast<LabelAbsPixelType> (abs(m_ForbiddenRegionLabel))) {
                    continue;
                }
                if (vLabelOfDefender != vLabelOfPropagatingRegion) {

                    ContourIndexType vNeighborIndex = vCurrentIndex;
                    vNeighborIndex += vOff;
                    
                    /// Tell the mother about the daughter:
                    (*aCandidateContainer)[vCurrentIndex].m_daughterIndices.push_back(vNeighborIndex);

                    CellListedHashMapIteratorType vContourParticleItr = aCandidateContainer->find(vNeighborIndex);
                    if (vContourParticleItr == aCandidateContainer->end()) {
                        /// create a new entry (a daughter), the contour point
                        /// has not been part of the contour so far.
                        ContourContainerValueType vOCCValue;
                        vOCCValue.m_candidateLabel = vLabelOfPropagatingRegion;
                        vOCCValue.m_label = vLabelOfDefender;
                        vOCCValue.m_data = this->GetInput()->GetPixel(vNeighborIndex);
                        vOCCValue.m_isDaughter = true;
                        vOCCValue.m_isMother = false;
                        vOCCValue.m_processed = false;
                        vOCCValue.m_referenceCount = 0; 

                        /// Tell the daughter about the mother:
                        vOCCValue.m_motherIndices.push_back(vCurrentIndex);
                        //                        vOCCValue.m_modifiedCounter = 0;

                        vOCCValue.m_energyDifference = NumericTraits<EnergyDifferenceType>::max();
                        vOCCValue.m_externalEnergyDifference = NumericTraits<EnergyDifferenceType>::max();

                        (*aCandidateContainer)[vNeighborIndex] = vOCCValue;
 
                    } else { /// the point is already part of the candidate list
                        vContourParticleItr->second.m_isDaughter = true;
                        /// Tell the daughter about the mother (label does not matter here!):
                        vContourParticleItr->second.m_motherIndices.push_back(vCurrentIndex);

                    }
                }
            }
        }
    }

    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    AddMothersToCandidateList(CellListedHashMapType* aCandidateContainer) {
        /// Add all the mother points - this is copying the inner contour list.
        /// (Things get easier afterwards if this is done in advance.)
        CellListedHashMapIteratorType vPointIterator = m_InnerContourContainer.begin();
        CellListedHashMapIteratorType vPointsEnd = m_InnerContourContainer.end();
        for (; vPointIterator != vPointsEnd; ++vPointIterator) {
            ContourContainerKeyType vCurrentIndex = vPointIterator->first;
            ContourContainerValueType vVal = vPointIterator->second;

            if (m_UseFastEvolution && vVal.m_modifiedCounter != 0 &&
                    (m_iteration_counter) % vVal.m_modifiedCounter != 0) {
                continue;
            }

            vVal.m_candidateLabel = 0;
            vVal.m_referenceCount = 0; // doesn't matter for the BG
            vVal.m_isMother = true;
            vVal.m_isDaughter = false;
            vVal.m_processed = false;
            vVal.m_modifiedCounter = vVal.m_modifiedCounter + 1;
            vPointIterator->second.m_modifiedCounter++;
            vVal.m_energyDifference = 0;
            vVal.m_externalEnergyDifference = 0;
            vVal.m_motherIndices.clear(); // this is necessary!
            vVal.m_daughterIndices.clear(); // this is necessary!!

            (*aCandidateContainer)[vCurrentIndex] = vVal;
        }
    }

    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    GetPropagatingCandidates(CellListedHashMapType * aCandidateContainer,
    std::vector<std::vector<int> >* aCoords,
    std::vector<LabelAbsPixelType>* aCandLabels,
    std::vector<unsigned int>* aReferenceCounts) {

        CellListedHashMapIteratorType vPointIterator = aCandidateContainer->begin();
        CellListedHashMapIteratorType vPointsEnd = aCandidateContainer->end();

        for (; vPointIterator != vPointsEnd; ++vPointIterator) {

            ContourContainerKeyType vCurrentIndex = vPointIterator->first;
            ContourContainerValueType vVal = vPointIterator->second;
            
            /// If this point is not a daughter, the energy has been calculated
            /// already (shrinking step).
            if (!vPointIterator->second.m_isDaughter) {
                continue;
            }

            std::set<LabelPixelType> vTestedLabels;

            /// Iterate the mothers
            typedef typename ContourParticleType::MotherIndexListType::iterator
            MotherListIteratorType;
            MotherListIteratorType vDMIt = vPointIterator->second.m_motherIndices.begin();
            MotherListIteratorType vDMItEnd = vPointIterator->second.m_motherIndices.end();
            for (; vDMIt != vDMItEnd; ++vDMIt) {
                ContourIndexType vDM = *vDMIt;
                
                ContourContainerValueType vMother = aCandidateContainer->find(vDM)->second;
                LabelAbsPixelType vPropagatingLabel = vMother.m_label;

                if (vTestedLabels.find(vPropagatingLabel) == vTestedLabels.end()) {
                    vTestedLabels.insert(vPropagatingLabel);

                    unsigned int vNbMothersOfPropLabel = 0;
                    MotherListIteratorType vMothersIt = vPointIterator->second.m_motherIndices.begin();
                    for (; vMothersIt != vDMItEnd; ++vMothersIt) {
                        
                        // Accessing the []-operator (in a map-iteration loop)
                        // might invalidate the iterator! 
                        // Don't use (*aCandidateContainer)[*vMothersIt];
                        ContourContainerValueType vMother2 = aCandidateContainer->find(*vMothersIt)->second;
                        if (vMother2.m_label == vPropagatingLabel) {
                            vNbMothersOfPropLabel++;
                        }
                    }
                    (*aReferenceCounts).push_back(vNbMothersOfPropLabel);

                    for (unsigned int vJ = 0; vJ < m_Dim; vJ++) {
                        (*aCoords)[vJ].push_back(vPointIterator->first[vJ]);
                    }
                    (*aCandLabels).push_back(vPropagatingLabel);
                }
            }
        }
    }

    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    MinimizeEnergy(
    CellListedHashMapType* aCandidateContainer,
    std::vector<std::vector<GPU_IntType> >* aCoords,
    std::vector<LabelAbsPixelType>* aCandLabels,
    std::vector<EnergyDifferenceType>* aInternalEnergy,
    std::vector<EnergyDifferenceType>* aExternalEnergy,
    std::vector<EnergyEventType>* aMerge,
    std::vector<unsigned int>* aReferenceCounts) {

        ContourIndexType vCIndex;
        unsigned int vNbParticles = (*aCoords)[0].size();

        typedef std::set<ContourIndexType> CIndexSetType;
        CIndexSetType vMergingSeeds;

        for (unsigned int vP = 0; vP < vNbParticles; vP++) { // particle
            for (unsigned int vD = 0; vD < m_Dim; vD++) {
                vCIndex[vD] = (*aCoords)[vD][vP];
            }
         
            ContourContainerValueType& vCurrentBestParticle = (*aCandidateContainer)[vCIndex];

            EnergyDifferenceType vIntEnergy = (*aInternalEnergy)[vP];
            EnergyDifferenceType vExtEnergy = (*aExternalEnergy)[vP];
            LabelAbsPixelType vCandidateLabel = (*aCandLabels)[vP];
            unsigned int vReferenceCount = (*aReferenceCounts)[vP];

            if(vIntEnergy != vIntEnergy || vExtEnergy != vExtEnergy) {
                assert(!"A energy value is is NaN.");
            }
                
            if (vExtEnergy + vIntEnergy <
                    vCurrentBestParticle.m_energyDifference + 
                    vCurrentBestParticle.m_externalEnergyDifference) {

                vCurrentBestParticle.m_externalEnergyDifference = vExtEnergy;
                vCurrentBestParticle.m_energyDifference = vIntEnergy;
                vCurrentBestParticle.m_candidateLabel = vCandidateLabel;
                vCurrentBestParticle.m_referenceCount = vReferenceCount;

                if ((*aMerge)[vP] && 
                        vCurrentBestParticle.m_label != 0 && 
                        vCurrentBestParticle.m_candidateLabel != 0) {
                    vMergingSeeds.insert(vCIndex);
                } else {
                    vMergingSeeds.erase(vCIndex); // might be expensive
                }
            }
            
        }

        typename CIndexSetType::iterator vSeedIt = vMergingSeeds.begin();
        typename CIndexSetType::iterator vSeedItEnd = vMergingSeeds.end();
        for (; vSeedIt != vSeedItEnd; ++vSeedIt) {
            /// add an entry in the global seed point list:
            ContourContainerValueType& vCurrentBestParticle = (*aCandidateContainer)[*vSeedIt];
            LabelAbsPixelType vL1 = vCurrentBestParticle.m_candidateLabel;
            LabelAbsPixelType vL2 = vCurrentBestParticle.m_label;

            CompetingRegionsPairType vPair;
            vPair[0] = (vL1 < vL2) ?
                    vL1 : vL2;
            vPair[1] = (vL1 > vL2) ?
                    vL1 : vL2;

            m_CompetingRegionsMap[vPair] = *vSeedIt;

            /// Ensure the point does not move since we'd like to merge
            /// here. Todo so, we set the energy to a large value. We are not 
            /// allowed to remove (erase) the point from the candidate list yet
            /// because the topology dependency processing would fail.
            (*aCandidateContainer)[*vSeedIt].m_energyDifference = NumericTraits<EnergyDifferenceType>::max();
        }

    }

    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    ComputeEnergiesForMothers(CellListedHashMapType* aCandidateContainer) {

        CellListedHashMapIteratorType vPointIterator = aCandidateContainer->begin();
        CellListedHashMapIteratorType vPointsEnd = aCandidateContainer->end();

        for (; vPointIterator != vPointsEnd; ++vPointIterator) {
            ContourContainerValueType& vParticleReference = vPointIterator->second;
            if (vParticleReference.m_isMother) {
                ContourContainerKeyType vCurrentIndex = vPointIterator->first;
                vParticleReference.m_energyDifference =
                        CalculateInternalEnergyDifferenceForLabel(vCurrentIndex, &(vPointIterator->second), 0);
                vParticleReference.m_externalEnergyDifference =
                        CalculateExternalEnergyDifferenceForLabel(vCurrentIndex, &(vPointIterator->second), 0).first;
                /// logically merge energy is missing here, but we treat this
                /// separately (not as a member of particle) and for daughters only.
            }
        }
    }

    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    ComputeEnergiesForDaughters(
    CellListedHashMapType* aCandidateContainer,
    std::vector<std::vector<GPU_IntType> >* aCoords,
    std::vector<LabelAbsPixelType>* aCandLabels,
    std::vector<EnergyDifferenceType>* aInternalEnergy,
    std::vector<EnergyDifferenceType>* aExternalEnergy,
    std::vector<EnergyEventType>* aMerge) {

        ContourIndexType vCurrentIndex;
        unsigned int vNbParticles = (*aCoords)[0].size();
        aInternalEnergy->resize(vNbParticles);
        aExternalEnergy->resize(vNbParticles);
        
        for (unsigned int vP = 0; vP < vNbParticles; vP++) { // particle
            for (unsigned int vD = 0; vD < m_Dim; vD++) { // dimension
                vCurrentIndex[vD] = (*aCoords)[vD][vP];
            }

            // Get a copy of the particle at the index of interest which is 
            // already in the candidate list
            ContourContainerValueType vParticle = (*aCandidateContainer)[vCurrentIndex];

            (*aInternalEnergy)[vP] = CalculateInternalEnergyDifferenceForLabel(
                    vCurrentIndex, &vParticle, (*aCandLabels)[vP]);

            EnergyDifferenceAndMergePairType vV = CalculateExternalEnergyDifferenceForLabel(
                    vCurrentIndex, &vParticle, (*aCandLabels)[vP]);
            (*aExternalEnergy)[vP] = vV.first;

            (*aMerge)[vP] = vV.second;
        }
    }


    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    RebuildAndOptimizeCandidateList(CellListedHashMapType* aCandidateContainer) {
        
        itk::TimeProbe vTimerRebuilding;
        vTimerRebuilding.Start();
        aCandidateContainer->clear();

        /// Copy the current contour, candidate labels are set to 0 (shrink)
        AddMothersToCandidateList(aCandidateContainer);

        /// Add daughter particles (candidate label is not important yet) in
        /// order to determine mother-daughter relations.
        AddDaughtersToCandidateList(aCandidateContainer);

        unsigned int vPreAllocSize = static_cast<unsigned int> (1.5f * 
        m_InnerContourContainer.size());

        std::vector<EnergyDifferenceType> vRMInternalEnergy;
        std::vector<EnergyDifferenceType> vRMExternalEnergy;
        std::vector<EnergyEventType> vRMMergeEnergy;
        std::vector<LabelAbsPixelType> vRMCandidateLabel;
        std::vector<unsigned int> vRMReferenceCounts;
        std::vector<std::vector<int> > vRMCoordinates(m_Dim);

        vRMCandidateLabel.reserve(vPreAllocSize);
        vRMReferenceCounts.reserve(vPreAllocSize);

        for (unsigned int i = 0; i < m_Dim; i++) {
            vRMCoordinates[i] = std::vector<int>();
            vRMCoordinates[i].reserve(vPreAllocSize);
        }

        /// Get the coordinates and candidate labels for all daughters (resp.
        /// for all propagating particles).
        GetPropagatingCandidates(aCandidateContainer, &vRMCoordinates, 
                &vRMCandidateLabel, &vRMReferenceCounts);


        vRMInternalEnergy.resize(vRMCandidateLabel.size());
        vRMExternalEnergy.resize(vRMCandidateLabel.size());
        vRMMergeEnergy.resize(vRMCandidateLabel.size());

        vTimerRebuilding.Stop();
        m_TimeHistRebuilding.push_back(vTimerRebuilding.GetTotal());


        /// Since the mother-daughter relationships are clear now, we can start
        /// to compute the energies. First, the differences in energy for the
        /// mother points when changing to the BG are calculated. The candidate
        /// labels were already set.
        itk::TimeProbe vTimerEnergyComputation;
        vTimerEnergyComputation.Start();

        if (m_GPUEnergy.IsNotNull()) {

            /// Define GPU arrays
            unsigned int **gInCoordinates;
            unsigned int* gInToLabel;
            float* gOutInternalEnergy;
            float* gOutExternalEnergy;
            unsigned int* gOutMerge;

            /// Allocate and Initialize GPU arrays // keep as is.
            unsigned int vPaddedNbParticles = InitializeGPUArrays(
                    aCandidateContainer, // mothers
                    &vRMCoordinates, // daughter coords
                    &vRMCandidateLabel, // daughter candidate labels
                    &gInCoordinates,
                    &gInToLabel,
                    &gOutInternalEnergy,
                    &gOutExternalEnergy,
                    &gOutMerge);

            m_GPUEnergy->EvaluateEnergyDifferences(
                    gInCoordinates,
                    gInToLabel,
                    gOutInternalEnergy,
                    gOutExternalEnergy,
                    gOutMerge,
                    vPaddedNbParticles);

            CopyGPUResults(
                    aCandidateContainer, // mothers
                    &vRMInternalEnergy,
                    &vRMExternalEnergy,
                    &vRMMergeEnergy,
                    gOutInternalEnergy,
                    gOutExternalEnergy,
                    gOutMerge);


            DistroyGPUArrays(
                    gInCoordinates,
                    gInToLabel,
                    gOutInternalEnergy,
                    gOutExternalEnergy,
                    gOutMerge);

        } else { // CPU based energies - standard procedure.

            ComputeEnergiesForMothers(aCandidateContainer); // only compute the mothers

            ComputeEnergiesForDaughters(
                    aCandidateContainer,
                    &vRMCoordinates,
                    &vRMCandidateLabel,
                    &vRMInternalEnergy,
                    &vRMExternalEnergy,
                    &vRMMergeEnergy);
        }

        vTimerEnergyComputation.Stop();
        m_TimeHistEnergyComputation.push_back(vTimerEnergyComputation.GetTotal());

        /// It remains to do the energy calculation for daughter points and the
        /// corresponding OPTIMIZATION (replace the candidate label for the label
        /// of lowest energy). We also increase the reference count of the mothers with the
        /// current best label.
        itk::TimeProbe vTimerMinimization;
        vTimerMinimization.Start();
        MinimizeEnergy(
                aCandidateContainer,
                &vRMCoordinates,
                &vRMCandidateLabel,
                &vRMInternalEnergy,
                &vRMExternalEnergy,
                &vRMMergeEnergy,
                &vRMReferenceCounts);
        vTimerMinimization.Stop();
        m_TimeHistMinimization.push_back(vTimerMinimization.GetTotal());
    }


    /// The function relabels touching regions. These 2 or more regions are
    /// fused into a new region. Both regions are at position of aIndex
    /// or adjacent to it.
    /// RelabelAdjacentRegionsAfterToplogicalChange expects the label image to
    /// be updated already: both methods are connected via the seedpoint. This
    /// method may be used to fuse 2 or more regions in the region competition mode.
    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    RelabelRegionsAfterFusions(LabelImagePointerType aLabelImage,
    InputImageIndexType aIndex, LabelAbsPixelType aL1) {
        typename MultipleThsFunctionType::Pointer vMultiThsFunction =
                MultipleThsFunctionType::New();
        vMultiThsFunction->SetInputImage(aLabelImage);

        /// if one of the labels is involved in another fusion, we'll also relabel
        /// this third/fourth... region
        typename CompetingRegionsMapType::iterator vCRit = m_CompetingRegionsMap.begin();
        typename CompetingRegionsMapType::iterator vCRend = m_CompetingRegionsMap.end();

        /// vLabelsToCheck contains all labels whose neighboring regions we still
        /// have to check.
        std::list<unsigned int> vLabelsToCheck;
        /// CheckedLabels contains all labels that are already in the multithsfunction.
        std::set<unsigned int> vCheckedLabels;


        /// Insert the first element (then iterate)
        vLabelsToCheck.push_back(aL1);
        vMultiThsFunction->AddThresholdBetween(aL1, aL1);
        vMultiThsFunction->AddThresholdBetween(-aL1, -aL1);
        vCheckedLabels.insert(aL1);

        std::list<unsigned int>::iterator vLabelsToCheckIt = vLabelsToCheck.begin();
        for (; vLabelsToCheckIt != vLabelsToCheck.end(); ++vLabelsToCheckIt) {
            vCRit = m_CompetingRegionsMap.begin();
            for (; vCRit != vCRend; ++vCRit) {
                LabelAbsPixelType vLabel1 = (vCRit->first)[0];
                LabelAbsPixelType vLabel2 = (vCRit->first)[1];

                if (vLabel1 == *vLabelsToCheckIt &&
                        vCheckedLabels.find(vLabel2) == vCheckedLabels.end()) {
                    vMultiThsFunction->AddThresholdBetween(vLabel2, vLabel2);
                    vMultiThsFunction->AddThresholdBetween(-vLabel2, -vLabel2);
                    vCheckedLabels.insert(vLabel2);
                    vLabelsToCheck.push_back(vLabel2);
                }
                if (vLabel2 == *vLabelsToCheckIt &&
                        vCheckedLabels.find(vLabel1) == vCheckedLabels.end()) {
                    vMultiThsFunction->AddThresholdBetween(vLabel1, vLabel1);
                    vMultiThsFunction->AddThresholdBetween(-vLabel1, -vLabel1);
                    vCheckedLabels.insert(vLabel1);
                    vLabelsToCheck.push_back(vLabel1);
                }
            }
        }

        if (vMultiThsFunction->EvaluateAtIndex(aIndex)) {
            ForestFire(aLabelImage, aIndex, vMultiThsFunction, m_MaxNLabels++);
        }
    }


    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    RelabelRegionsAfterSplit(
    LabelImagePointerType aLabelImage,
    InputImageIndexType aIndex,
    unsigned int aLabel) {

        /// Maybe it has been relabelled already. It is not necessary but should
        /// speed up
        if (abs(m_LabelImage->GetPixel(aIndex)) == aLabel) {

            typename MultipleThsFunctionType::Pointer vMultiThsFunction = MultipleThsFunctionType::New();
            vMultiThsFunction->SetInputImage(aLabelImage);
            vMultiThsFunction->AddThresholdBetween(aLabel, aLabel);
            vMultiThsFunction->AddThresholdBetween(-aLabel, -aLabel);

            ForestFire(m_LabelImage, aIndex, vMultiThsFunction, m_MaxNLabels++);
        }
    }

    /// The function relabels a region starting from the position aIndex.
    /// This method assumes the label image to be updated. It is used to relabel
    /// a region that was split by another region (maybe BG region).
    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    RegisterSeedsAfterSplit(
    LabelImagePointerType aLabelImage,
    ContourIndexType aIndex,
    LabelAbsPixelType aLabel,
    CellListedHashMapType* aCandidateContainer) {

        for (unsigned int vI = 0; vI < m_NeighborhoodSize_FG_Connectivity; vI++) {

            InputImageOffsetType vOff = m_NeighborsOffsets_FG_Connectivity[vI];
            LabelPixelType vLabel = aLabelImage->GetPixel(aIndex + vOff);

            if (abs(vLabel) == aLabel) {
                SeedType vSeed;
                for (unsigned int vD = 0; vD < m_Dim; vD++) {
                    vSeed.first[vD] = aIndex[vD] + vOff[vD];
                }
                vSeed.second = aLabel;
                m_Seeds.insert(vSeed);


                /// At the position where we put the seed, inform the particle
                /// that it has to inform its neighbor in case it moves (if there
                /// is a particle at all at this spot; else we don't have a problem
                /// because the label will not move at the spot and therefore the
                /// seed will be effective).
                CellListedHashMapIteratorType vCorrespParticle = aCandidateContainer->find(vSeed.first);
                if (vCorrespParticle != aCandidateContainer->end())
                    vCorrespParticle->second.m_processed = true;
            }
        }
    }

    /*
     * RelabelAllAdjacentRegionsAfterTopologicalChange assumes the label image
     * to be updated: The non-simple point pixel already has the new label A.
     * When splitting (A=0), each adjacent region gets a new label and new stats.
     * When merging, all regions adjacent regions get all the same label and new
     * stats.
     *
     * Thus, for competing regions, if they should not merge, this method
     * must not be invoked.
     * 
     * The method is generally not appropriate for RC mode.
     *
     * When splitting, all regions (except the propagating region) adjacent to
     * the non-simple point get relabeled. To avoid several passes, the method
     * only relabels regions that are smaller than aMaxLabel.
     */
    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    RelabelAllAdjacentRegionsAfterTopologicalChange(
    LabelImagePointerType aLabelImage,
    InputImageIndexType aIndex,
    LabelAbsPixelType aMaxLabel) {
        typename MultipleThsFunctionType::Pointer vMultiThsFunction = MultipleThsFunctionType::New();
        vMultiThsFunction->SetInputImage(aLabelImage);
        for (unsigned int vI = 0; vI < m_NeighborhoodSize_FG_Connectivity; vI++) {
            InputImageOffsetType vOff = m_NeighborsOffsets_FG_Connectivity[vI];
            LabelAbsPixelType vLabelAbs = static_cast<unsigned int> (
                    abs(aLabelImage->GetPixel(aIndex + vOff)));
            if (vLabelAbs != 0 &&
                    vLabelAbs != static_cast<unsigned int> (m_ForbiddenRegionLabel) &&
                    vLabelAbs < aMaxLabel) {
                vMultiThsFunction->AddThresholdBetween(vLabelAbs, vLabelAbs);
                vMultiThsFunction->AddThresholdBetween(-vLabelAbs, -vLabelAbs);
            }
        }
        ForestFire(aLabelImage, aIndex, vMultiThsFunction);
    }

    /**
     * The method ForestFire creates new regions for all regions adjacent to
     * aIndex AND are encoded in the argument MultiThsFunction.
     * The method assumes the voxel causing the topological change to be updated
     * already. The method can be used for merging and splitting regions.
     */
    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    ForestFire(
    LabelImagePointerType aLabelImage,
    InputImageIndexType aIndex,
    MultipleThsFunctionPointerType aMultiThsFunctionPtr,
    LabelAbsPixelType aNewLabel) {

        typedef std::set<LabelAbsPixelType> LabelAbsPxSetType;
        LabelAbsPxSetType vVisitedOldLabels;

        typedef FloodFilledImageFunctionConditionalIterator<
                LabelImageType, MultipleThsFunctionType> LabelIteratorType;
        LabelIteratorType vLit = LabelIteratorType(aLabelImage, aMultiThsFunctionPtr, aIndex);

        typedef std::set<ContourIndexType> ContourIndexSetType;
        ContourIndexSetType vSetOfAncientContourIndices;
        
        for (vLit.GoToBegin(); !vLit.IsAtEnd(); ++vLit) {
            LabelPixelType vLabelValue = vLit.Get();
            InputImageIndexType vCurrentIndex = vLit.GetIndex();
            InputPixelType vImageValue = this->GetInput()->GetPixel(vCurrentIndex);

            // the visited labels statistics will be removed later.
            vVisitedOldLabels.insert(abs(vLabelValue));


            vLit.Set(aNewLabel);
            if (vLabelValue < 0) {
                // TODO: GetIndex() in the floodFillIterator seems to be very expensive.
                //       Maybe there is a more convenient way: When ContourContainer will
                //       be semistructured, then just iterate through the object boundaries.

                // TODO: Cast from index to contourIndex doesn't work.
                ContourIndexType vCurrentCIndex;
                for (unsigned int vD = 0; vD < m_Dim; vD++) {
                    vCurrentCIndex[vD] = vCurrentIndex[vD];
                }

                vSetOfAncientContourIndices.insert(vCurrentCIndex);
            }

            for(unsigned int vI = 0; vI < m_ExtEnergies.size(); vI++) {
                m_ExtEnergies[vI]->AddPoint(vCurrentIndex, aNewLabel, vImageValue);
            }
            for(unsigned int vI = 0; vI < m_IntEnergies.size(); vI++) {
                m_IntEnergies[vI]->AddPoint(vCurrentIndex, aNewLabel);
            }
            for(unsigned int vI = 0; vI < m_PIEnergies.size(); vI++) {
                m_PIEnergies[vI]->AddPoint(vCurrentIndex, aNewLabel);
            }
            if(m_GPUEnergy.IsNotNull()) {
                m_GPUEnergy->AddPoint(vCurrentIndex, aNewLabel, vImageValue);
            }
        }


        /// Delete the contour points that are not needed anymore:
        typename ContourIndexSetType::iterator vCPIt = vSetOfAncientContourIndices.begin();
        typename ContourIndexSetType::iterator vCPItEnd = vSetOfAncientContourIndices.end();
        for (; vCPIt != vCPItEnd; ++vCPIt) {
            ContourIndexType vCurrentCIndex = *vCPIt;

            if (IsBoundaryPoint(vCurrentCIndex)) {
                ContourContainerValueType& vPoint = m_InnerContourContainer[vCurrentCIndex];
                vPoint.m_label = aNewLabel;
                vPoint.m_modifiedCounter = 0;

                m_LabelImage->SetPixel(vCurrentCIndex, -aNewLabel);
            } else {
                m_InnerContourContainer.erase(vCurrentCIndex);
            }
        }

        /// Clean up the statistics of non valid regions.
        LabelAbsPxSetType::iterator vVisitedIt = vVisitedOldLabels.begin();
        LabelAbsPxSetType::iterator vVisitedItEnd = vVisitedOldLabels.end();
        for (; vVisitedIt != vVisitedItEnd; ++vVisitedIt) {
            FreeLabelStatistics(*vVisitedIt);
        }
    }

    
    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    FreeLabelStatistics(LabelAbsPixelType aLabelAbs) {

        for(unsigned int vI = 0; vI < m_ExtEnergies.size(); vI++) {
            m_ExtEnergies[vI]->KillRegion(aLabelAbs);
        }
        for(unsigned int vI = 0; vI < m_IntEnergies.size(); vI++) {
            m_IntEnergies[vI]->KillRegion(aLabelAbs);
        }
        for(unsigned int vI = 0; vI < m_PIEnergies.size(); vI++) {
            m_PIEnergies[vI]->KillRegion(aLabelAbs);
        }
        if(m_GPUEnergy.IsNotNull()) {
            m_GPUEnergy->KillRegion(aLabelAbs);
        }
                
    }


    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    PrepareEnergyCaluclationForEachIteration() {

        for(unsigned int vI = 0; vI < m_ExtEnergies.size(); vI++) {
            m_ExtEnergies[vI]->PrepareEnergyCalculationForIteration();
        }
        for(unsigned int vI = 0; vI < m_IntEnergies.size(); vI++) {
            m_IntEnergies[vI]->PrepareEnergyCalculationForIteration();
        }
        for(unsigned int vI = 0; vI < m_PIEnergies.size(); vI++) {
            m_PIEnergies[vI]->PrepareEnergyCalculationForIteration();
        }
        if(m_GPUEnergy.IsNotNull()) {
            m_GPUEnergy->PrepareEnergyCalculationForIteration();
        }
        
    }

    /**
     * Method is used to prepare special energy functions. It is called only once
     * in the beginning of the filter. Use PrepareEnergyCalculationForIteration()
     * to have a function that is called at every iteration.
     * For example for the deconvolution energy:
     * This method calculates the point spread function and creates the
     * 'ideal image'.
     * The method assumes that statistics for the regions have already been
     * calculated.
     */
    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    PrepareEnergyCaluclation() {
        
        for(unsigned int vI = 0; vI < m_ExtEnergies.size(); vI++) {
            m_ExtEnergies[vI]->PrepareEnergyCalculation();
        }
        for(unsigned int vI = 0; vI < m_IntEnergies.size(); vI++) {
            m_IntEnergies[vI]->PrepareEnergyCalculation();        
        }
        for(unsigned int vI = 0; vI < m_PIEnergies.size(); vI++) {
            m_PIEnergies[vI]->PrepareEnergyCalculation();        
        }
        if(m_GPUEnergy.IsNotNull()) {
            m_GPUEnergy->PrepareEnergyCalculation();
        }
    }


    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    CopyRegionsAndAllocateOutput() {
        //        this->GetOutput()->SetRequestedRegion(this->GetInput()->GetRequestedRegion());
        this->GetOutput()->SetBufferedRegion(this->GetInput()->GetBufferedRegion());
        //        this->GetOutput()->SetLargestPossibleRegion(this->GetInput()->GetLargestPossibleRegion());

        this->GetOutput()->Allocate();
    }

    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    InitializeLabelImageAndContourContainer(InitImagePointerType aInitImage) {
        /**
         * Create the labeling function and connect it to the output of <code>this</code>.
         */

        int vMins[m_Dim];
        int vMaxs[m_Dim];
        for (unsigned int vD = 0; vD < m_Dim; vD++) {
            vMins[vD] = this->GetInput()->GetBufferedRegion().GetIndex()[vD];
            vMaxs[vD] = vMins[vD] + this->GetInput()->GetBufferedRegion().GetSize()[vD] - 1;
        }

        typename InitImageType::Pointer vForbiddenRegionImage;

        if (m_UseForbiddenRegion) {
            vForbiddenRegionImage = this->GetForbiddenRegionInput();
        } else {
            vForbiddenRegionImage = InitImageType::New();
            vForbiddenRegionImage->SetRegions(aInitImage->GetBufferedRegion());
            vForbiddenRegionImage->Allocate();
            vForbiddenRegionImage->FillBuffer(0);
            vForbiddenRegionImage->SetSpacing(this->GetDataInput()->GetSpacing());
        }

        typedef ImageRegionIteratorWithIndex<InitImageType> ForbImageRegItWithIndexType;
        ForbImageRegItWithIndexType vForbIt(vForbiddenRegionImage,
                vForbiddenRegionImage->GetBufferedRegion());
        for (vForbIt.GoToBegin(); !vForbIt.IsAtEnd(); ++vForbIt) {
            InitImageIndexType vIndex = vForbIt.GetIndex();
            for (unsigned int vD = 0; vD < m_Dim; vD++) {
                if (vIndex[vD] == vMins[vD] || vIndex[vD] == vMaxs[vD]) {
                    vForbIt.Set(1);
                    continue;
                }
            }
        }

        typedef itk::BinaryThresholdImageFilter<InitImageType, LabelImageType>
                BinaryThresholdImageFilterType;
        typename BinaryThresholdImageFilterType::Pointer vBinaryThresholdFilter =
                BinaryThresholdImageFilterType::New();

        vBinaryThresholdFilter->SetLowerThreshold(0);
        vBinaryThresholdFilter->SetUpperThreshold(0);
        vBinaryThresholdFilter->SetInsideValue(1);
        vBinaryThresholdFilter->SetOutsideValue(0);
        vBinaryThresholdFilter->SetInput(vForbiddenRegionImage);


        typedef itk::MaskImageFilter<InitImageType, LabelImageType, LabelImageType>
                MaskImageFilterType;
        typename MaskImageFilterType::Pointer vMaskFilter = MaskImageFilterType::New();

        vMaskFilter->SetOutsideValue(0);
        vMaskFilter->SetInput1(aInitImage);
        vMaskFilter->SetInput2(vBinaryThresholdFilter->GetOutput());

        typedef itk::ScalarConnectedComponentImageFilter<LabelImageType, LabelImageType>
                ScalarConnectedCompFilterType;
        typename ScalarConnectedCompFilterType::Pointer vScalarConnCompFilter =
                ScalarConnectedCompFilterType::New();
        vScalarConnCompFilter->SetInput(vMaskFilter->GetOutput());
        vScalarConnCompFilter->SetDistanceThreshold(0);
        vScalarConnCompFilter->SetMaskImage(vMaskFilter->GetOutput());

        // TODO: what is fully connected on: itk doc is not clear!
        //        vScalarConnCompFilter->FullyConnectedOn();

        vScalarConnCompFilter->Update();
        m_LabelImage = vScalarConnCompFilter->GetOutput();


        
        /*
         * Get the contour of the labels
         */
        typedef itk::LabelContourImageFilter<LabelImageType, LabelImageType> LabelContourFilterType;
        typename LabelContourFilterType::Pointer vLabelContourFilter = LabelContourFilterType::New();
        vLabelContourFilter->SetBackgroundValue(0);
        vLabelContourFilter->SetInput(vScalarConnCompFilter->GetOutput());
        vLabelContourFilter->SetInput(m_LabelImage);

        vLabelContourFilter->Update();

        /*
         * Set up the containers for the label and the input image
         */
        typedef itk::ImageRegionConstIterator<InputImageType> DataImageIteratorType;
        DataImageIteratorType vDataImageIterator(
                this->GetDataInput(), this->GetInput()->GetBufferedRegion());

        /*
         * Find the cell size for the cell-listed hashmap.
         * Clear and (re-)fill the containers that store the contour/front.
         */
        m_InnerContourContainer.clear();
        unsigned int vCellSize[m_Dim];
        unsigned int vRegionSize[m_Dim];
        int vRegionIndex[m_Dim];
                       
        for (unsigned int vD = 0; vD < m_Dim; vD++) {
            vRegionSize[vD] = static_cast<unsigned int> (abs(vMaxs[vD] - vMins[vD]));
            vRegionIndex[vD] = vMins[vD];
            
            // set the minimum cell size
            vCellSize[vD] = m_MinimalCellSize[vD];
            
            // check if we need to enlarge the cell size (requested by a energy)
            for(unsigned vPI = 0; vPI < m_PIEnergies.size(); vPI++) {
                if(vCellSize[vD] < m_PIEnergies[vPI]->GetMinimalIteratorSize()[vD]) {
                    vCellSize[vD] = m_PIEnergies[vPI]->GetMinimalIteratorSize()[vD];
                }
            }
            
            // check if we need to enlarge the cell size because of the 
            // external energy smoothing (Sobolev stuff):
            if (m_UseSobolevGradients) {
                if (3 * m_SobolevKernelSigma[vD] > vCellSize[vD]) {
                    vCellSize[vD] = 3 * (m_SobolevKernelSigma[vD]);
                }
            }
        }
        m_InnerContourContainer.SetCellAndRegionSize(vCellSize, vRegionSize, vRegionIndex);
        m_Candidates.SetCellAndRegionSize(vCellSize, vRegionSize, vRegionIndex);
        
        /**
         * Fix the boundary and the forbidden region
         */
        typedef itk::ImageRegionIteratorWithIndex<LabelImageType> LabelImageIteratorType;
        LabelImageIteratorType vLabelIt(m_LabelImage, m_LabelImage->GetBufferedRegion());
        if (m_UseForbiddenRegion) {
            /// Put the forbidden region into the label image:
            typedef itk::ImageRegionConstIterator<InitImageType> InitImageConstIteratorType;
            InitImageConstIteratorType vForbIt(
                    this->GetForbiddenRegionInput(), this->GetForbiddenRegionInput()->GetBufferedRegion());

            for (vForbIt.GoToBegin(), vLabelIt.GoToBegin();
                    !vLabelIt.IsAtEnd();
                    ++vLabelIt, ++vForbIt) {
                if (vForbIt.Get() != NumericTraits<InitPixelType>::Zero) {
                    vLabelIt.Set(m_ForbiddenRegionLabel);
                }
            }
        }
        
        for(vLabelIt.GoToBegin(); !vLabelIt.IsAtEnd(); ++vLabelIt) {
            LabelImageIndexType vIndex = vLabelIt.GetIndex();
                        bool vBoundaryPx = false;
            for (unsigned int vD = 0; vD < m_Dim; vD++) {
                if (vIndex[vD] == vMins[vD] || vIndex[vD] == vMaxs[vD]) {
                    vBoundaryPx = true;
                }
            }

            if (vBoundaryPx) {
                m_LabelImage->SetPixel(vIndex, m_ForbiddenRegionLabel);
                continue;
            }
        }
        
        /// Insert the index to the contour container, set the label image to a
        /// negative value and initialize the memory image:
        LabelImageIteratorType vLabelContourImageIterator(
        vLabelContourFilter->GetOutput(),
                vLabelContourFilter->GetOutput()->GetBufferedRegion());
        for (vLabelContourImageIterator.GoToBegin(), vDataImageIterator.GoToBegin();
                !vLabelContourImageIterator.IsAtEnd();
                ++vLabelContourImageIterator, ++vDataImageIterator) {
            
            ContourContainerKeyType vIndex;
            //TODO: cast somehow instead of copying!??
            for (unsigned int i = 0; i < InputImageType::ImageDimension; i++) {
                vIndex[i] = vLabelContourImageIterator.GetIndex()[i];
            }
            
            if (vLabelContourImageIterator.Get() != 0 &&
                    vLabelContourImageIterator.Get() != m_ForbiddenRegionLabel) {
                LabelAbsPixelType vAbsLabel = abs(m_LabelImage->GetPixel(vIndex));
                
                /// Fill the container that stores the inner contour: these pixel still
                /// belong to the object.
                ContourContainerValueType vVal;
                vVal.m_label = vAbsLabel;
                vVal.m_data = vDataImageIterator.Get();
                m_InnerContourContainer[vIndex] = vVal;
                
                /// To distinguish contour-pixel from pixel inside the object, the contour
                /// pixel have negative labels in the labelimage.
                m_LabelImage->SetPixel(vIndex, -vAbsLabel);
                m_MemoryImage->SetPixel(vIndex, m_MemoryLength);
            }
        }
        
    }
    
    /**
     * The method sets up the statistics of energies. It 
     * iterates the label and data image and reports each point to each energy
     * registered.
     */
    template <class TInputImage, class TInitImage, class TOutputImage >
    void FrontsCompetitionImageFilter<TInputImage, TInitImage, TOutputImage>::
    RenewStatistics(LabelImagePointerType aInitImage) {

        ImageRegionConstIteratorWithIndex<InputImageType> vDataIt(this->GetDataInput(),
                this->GetDataInput()->GetBufferedRegion());
        ImageRegionConstIterator<LabelImageType> vLabelIt(m_LabelImage,
                m_LabelImage->GetBufferedRegion());
        m_MaxNLabels = 0;
        for (vDataIt.GoToBegin(), vLabelIt.GoToBegin();
                !vLabelIt.IsAtEnd(); ++vDataIt, ++vLabelIt) {
            LabelAbsPixelType vAbsLabel = static_cast<LabelAbsPixelType> (abs(vLabelIt.Get()));
            InputPixelType vData = vDataIt.Get();
            if (vAbsLabel != static_cast<LabelAbsPixelType> (m_ForbiddenRegionLabel)) {
                
                for(unsigned int vI = 0; vI < m_ExtEnergies.size(); vI++) {
                    m_ExtEnergies[vI]->AddPoint(vDataIt.GetIndex(), vAbsLabel, vData);
                }
                for(unsigned int vI = 0; vI < m_IntEnergies.size(); vI++) {
                    m_IntEnergies[vI]->AddPoint(vDataIt.GetIndex(), vAbsLabel);
                }
                for(unsigned int vI = 0; vI < m_PIEnergies.size(); vI++) {
                    m_PIEnergies[vI]->AddPoint(vDataIt.GetIndex(), vAbsLabel);
                }
                if(m_GPUEnergy.IsNotNull()) {
                    m_GPUEnergy->AddPoint(vDataIt.GetIndex(), vAbsLabel, vData);
                }
                
                if (static_cast<LabelAbsPixelType> (abs(vLabelIt.Get())) > m_MaxNLabels) {
                    m_MaxNLabels = static_cast<LabelAbsPixelType> (abs(vLabelIt.Get()));
                }
            }
        }

        m_MaxNLabels++; // this number points to the a free label.
    }



} //end namespace itk

#endif //__ITKFRONTSCOMPETITIONIMAGEFILTER_CXX_
    