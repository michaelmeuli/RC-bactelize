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

#include "itkRCGPUPiecewiseConstantSquareDistCurvatureRegEnergy.h"

namespace itk
{
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    RCGPUPiecewiseConstantSquareDistCurvatureRegEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::RCGPUPiecewiseConstantSquareDistCurvatureRegEnergy() {
        m_RegionMergingThreshold = 1;
        m_DataCoeff = 1;
        m_LengthCoeff = 0.04;
        m_BalloonCoeff = 0.01;
        m_ConstantOutwardFlowCoeff = 0;
        m_CurvatureMaskRadius = 4;
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void RCGPUPiecewiseConstantSquareDistCurvatureRegEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::PrintSelf(std::ostream & os, Indent indent) const
    {
        Superclass::PrintSelf(os, indent); 
        os << indent << "Region merging threshold: " << m_RegionMergingThreshold << std::endl;
        os << indent << "External term coefficient: " << m_DataCoeff << std::endl;
        os << indent << "Internal length term coefficient: " << m_LengthCoeff << std::endl;
        os << indent << "Balloon flow coefficient: " << m_BalloonCoeff << std::endl;
        os << indent << "Constant outward flow coeff: " << m_ConstantOutwardFlowCoeff << std::endl;
        os << indent << "Curvature estimation mask radius " << m_CurvatureMaskRadius << std::endl;
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void RCGPUPiecewiseConstantSquareDistCurvatureRegEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::PrepareEnergyCalculation(){
//        CheckGPUusability();

        if(this->m_DataImage.IsNull() || this->m_LabelImage.IsNull()) {
            itkExceptionMacro("Data or label image not defined for GPU energy.");
        }
        
        InitOpenCL();

        InitGPUMembers();

        PadDataImageAndCopyToGPUArray();
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void RCGPUPiecewiseConstantSquareDistCurvatureRegEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::PrepareEnergyCalculationForIteration(){

        PadLabelImageAndCopyToGPUArray();
        
    }
        
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void RCGPUPiecewiseConstantSquareDistCurvatureRegEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::InitGPUMembers() {
        for (unsigned int vD = 0; vD < LabelImageType::ImageDimension; vD++) {
            float vScaleDim = this->m_DataImage->GetSpacing()[0] /
            this->m_DataImage->GetSpacing()[vD];
            m_gCurvatureRadius[vD] = m_CurvatureMaskRadius * vScaleDim;
            m_gRegionRadius[vD] = m_PiecewiseSmoothMaskRadius * vScaleDim;
        }
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void RCGPUPiecewiseConstantSquareDistCurvatureRegEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::InitOpenCL() {
        itk_init_clPlatform(&m_gPlatform);
        itk_init_clDevice(m_gPlatform, &m_gDevice);

        itk_init_clContextAndQueue(m_gPlatform, m_gDevice, &m_gContext, &m_gCmd_queue);

        itk_create_clKernel(m_gContext, m_gKernels, m_gProgram, m_gDevice, 
                "gpu_kernels/internal_curvature.cl", "energy", LabelImageType::ImageDimension, 0);
        itk_create_clKernel(m_gContext, m_gKernels, m_gProgram, m_gDevice, 
                "gpu_kernels/external_ps_gauss.cl", "energy", LabelImageType::ImageDimension, 1);
    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void RCGPUPiecewiseConstantSquareDistCurvatureRegEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::PadLabelImageAndCopyToGPUArray() {

        typedef ConstantPadImageFilter<LabelImageType, LabelImageType> PadLabelImageFilterType;
        typename PadLabelImageFilterType::Pointer vPadLabelImageFilter = PadLabelImageFilterType::New();
        // pad with the forbidden label of the optimizer (todo:get it somehow).
        vPadLabelImageFilter->SetConstant(NumericTraits<LabelPixelType>::max()); 

        /// get the max of the curvature radius and the data-patch radius
        for (unsigned int vD = 0; vD < LabelImageType::ImageDimension; vD++) {
            m_gLabelPaddingSize[vD] = (m_gCurvatureRadius[vD] > m_gRegionRadius[vD]) ?
                    ceil(m_gCurvatureRadius[vD]) : ceil(m_gRegionRadius[vD]);
        }

        vPadLabelImageFilter->SetPadLowerBound(m_gLabelPaddingSize);
        vPadLabelImageFilter->SetPadUpperBound(m_gLabelPaddingSize);
        vPadLabelImageFilter->SetInput(this->m_LabelImage);
        vPadLabelImageFilter->Update();

        ImageRegionConstIterator<LabelImageType> vLabelIt(
                vPadLabelImageFilter->GetOutput(), vPadLabelImageFilter->GetOutput()->GetLargestPossibleRegion());
        vLabelIt.GoToBegin();
        int vII = 0;
        for (vII = 0; !vLabelIt.IsAtEnd(); ++vLabelIt, vII++) {
            m_gLabelImage[vII] = static_cast<GPU_IntType> (vLabelIt.Get());
        }

    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void RCGPUPiecewiseConstantSquareDistCurvatureRegEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::PadDataImageAndCopyToGPUArray() {

        typedef ConstantPadImageFilter<DataImageType, DataImageType> PadDataImageFilterType;
        typename PadDataImageFilterType::Pointer vPadDataImageFilter = PadDataImageFilterType::New();
        //vPadDataImageFilter->SetConstant(m_ForbiddenRegionLabel);
        vPadDataImageFilter->SetConstant(0.0);

        /// get the max
        for (unsigned int vD = 0; vD < LabelImageType::ImageDimension; vD++) {
            m_gDataPaddingSize[vD] = (m_gCurvatureRadius[vD] > m_gRegionRadius[vD]) ?
                    ceil(m_gCurvatureRadius[vD]) : ceil(m_gRegionRadius[vD]);
        }

        vPadDataImageFilter->SetPadLowerBound(m_gDataPaddingSize);
        vPadDataImageFilter->SetPadUpperBound(m_gDataPaddingSize);
        vPadDataImageFilter->SetInput(this->m_DataImage);
        vPadDataImageFilter->Update();

        unsigned int vS = 1;
        SizeType vPaddedSize =
                vPadDataImageFilter->GetOutput()->GetLargestPossibleRegion().GetSize();
        for (unsigned int vD = 0; vD < LabelImageType::ImageDimension; vD++) {
            vS *= vPaddedSize[vD];
        }

        m_gDataImage = new float[vS];
        m_gLabelImage = new GPU_IntType[vS];

        ImageRegionConstIterator<DataImageType> vIt(
                vPadDataImageFilter->GetOutput(), 
                vPadDataImageFilter->GetOutput()->GetLargestPossibleRegion());
        vIt.GoToBegin();
        for (int vI = 0; !vIt.IsAtEnd(); ++vIt, vI++) {
            m_gDataImage[vI] = static_cast<float> (vIt.Get());
        }
    }
    
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void 
    RCGPUPiecewiseConstantSquareDistCurvatureRegEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::EvaluateEnergyDifferences( //ComputeEnergiesForAllParticlesOnGPU(
    unsigned int **aGPUCoordinates,
    unsigned int* aGPUToLabel,
    float* aGPUOutInternalEnergy,
    float* aGPUOutExternalEnergy,
    unsigned int* aGPUOutMerge,
    unsigned int aPaddedNbCandidates) {

        SizeType imageSize = 
                this->m_DataImage->GetLargestPossibleRegion().GetSize();

        /// we just work with 3D. In case we're in 2D, the 3rd is just set to 1
        /// (which chorresponds to 1 slice).
        unsigned int* vGPUInTotGlobSize = new unsigned int[3];
        unsigned int* vGPUInOffset = new unsigned int[3];
        vGPUInTotGlobSize[2] = 1;
        vGPUInOffset[2] = 1;
        for (int i = 0; i < LabelImageType::ImageDimension; i++) {
            vGPUInOffset[i] = ((m_gLabelPaddingSize[i] > m_gDataPaddingSize[i]) ?
                    m_gLabelPaddingSize[i] : m_gDataPaddingSize[i]);
            vGPUInTotGlobSize[i] = imageSize[i] + 2 * vGPUInOffset[i];
        }

        unsigned int vGPUInGlobalSize1D = aPaddedNbCandidates; // + (512 - (aPaddedNbCandidates & 511));
        unsigned int vGPUInLocalSize1D = m_gWorkGroupSize;

        float vGPUInVolume = this->CalculateScaledSphereVolume(m_CurvatureMaskRadius);

        if (m_LengthCoeff != 0) {
            calculateInternalEnergyWithGPU(
                    m_gLabelImage, vGPUInTotGlobSize, m_gCurvatureRadius,
                    vGPUInOffset, vGPUInVolume, m_LengthCoeff,
                    aGPUOutInternalEnergy, aGPUOutMerge, aGPUCoordinates,
                    aGPUToLabel, vGPUInGlobalSize1D, vGPUInLocalSize1D);
        }

        int numOfLabels = m_Sums.size();//sizeof (m_Means) / sizeof (m_Means[0]); /// WHuuuaaa
        float* m_gMeans = (float*) calloc(numOfLabels, sizeof (float));
        float* m_gVars = (float*) calloc(numOfLabels, sizeof (float));
        for (int i = 0; i < numOfLabels; i++) {
            m_gMeans[i] = static_cast<float> (m_Sums[i]/m_Count[i]);
            m_gVars[i] = static_cast<float> (this->CalculateVariance(m_Sums_2[i],m_gMeans[i],m_Count[i]));
        }


        if (m_DataCoeff != 0) {
            calculateExternalEnergyWithGPU(
                    m_gLabelImage, m_gDataImage, vGPUInTotGlobSize,
                    m_DataCoeff, m_RegionMergingThreshold,
                    m_BalloonCoeff, m_ConstantOutwardFlowCoeff,
                    m_gRegionRadius, vGPUInOffset,
                    aGPUOutExternalEnergy, aGPUOutMerge, aGPUCoordinates,
                    aGPUToLabel, vGPUInGlobalSize1D, vGPUInLocalSize1D,
                    m_gMeans, m_gVars, numOfLabels);
        }

        free(m_gMeans);
        free(m_gVars);
        delete [] vGPUInTotGlobSize;
	delete [] vGPUInOffset;
    }

    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void 
    RCGPUPiecewiseConstantSquareDistCurvatureRegEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::calculateInternalEnergyWithGPU(
    GPU_IntType* aGPUInLabelImg, unsigned int* aGPUInTotGlobSize,
    float* aGPUInMaskRadius, unsigned int* aGPUInOffset,
    float aGPUInVolume, float aGPUInEnergyContourLengthCoeff,
    float* aGPUOutCurvatureFlow, unsigned int* aGPUOutTesting,
    unsigned int** aGPUInCandidate, unsigned int* aGPUInToLabel,
    unsigned int aGPUInGlobalSize1D, unsigned int aGPUInLocalSize1D) {


        cl_int err;

        unsigned int slice = aGPUInTotGlobSize[0] * aGPUInTotGlobSize[1];

        size_t buffer_size_1D = sizeof (GPU_UIntType) * aGPUInGlobalSize1D; // number of particles padded to a multiple of 512
        size_t local_buffer_size_1D = sizeof (GPU_UIntType) * aGPUInLocalSize1D; // workgroup size of uints
        size_t float_buffer_size_1D = sizeof (float) * aGPUInGlobalSize1D;

        // the memory objects (reference to the global memory objects):
        cl_mem labels_mem, candidate_x_mem, candidate_y_mem, candidate_z_mem,
                curvature_flow_mem, testing_mem, to_labels_mem;

        cl_image_format format; // characterizes the image in texture memory.
        format.image_channel_data_type = CL_SIGNED_INT32;
        format.image_channel_order = CL_R;

        // start to create the objects:
        cl_int imageerr;

        // using createImage will create an image in TEXTURE memory
        if (LabelImageType::ImageDimension == 2) {
            labels_mem = clCreateImage2D(m_gContext,
                    CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                    &format,
                    aGPUInTotGlobSize[0],
                    aGPUInTotGlobSize[1],
                    aGPUInTotGlobSize[0] * sizeof (int), // !! provide a row size
                    aGPUInLabelImg, // c++ array
                    &imageerr);
        } else if (LabelImageType::ImageDimension == 3) {
            labels_mem = clCreateImage3D(m_gContext,
                    CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                    &format,
                    aGPUInTotGlobSize[0],
                    aGPUInTotGlobSize[1],
                    aGPUInTotGlobSize[2],
                    aGPUInTotGlobSize[0] * sizeof (int),
                    slice * sizeof (int),
                    aGPUInLabelImg,
                    &imageerr);
        }

        // the following creates nice output in case we have an error
        itk_handle_clEnqueueCreateImage(imageerr);
        assert(imageerr == CL_SUCCESS);


        // Allocating buffer object (standard arrays; it will be in standard memory).
        //        masks_mem = clCreateBuffer(m_gContext, CL_MEM_READ_ONLY, mask_buffer_size, NULL, NULL);

        candidate_x_mem = clCreateBuffer(m_gContext, CL_MEM_READ_ONLY, buffer_size_1D, NULL, NULL);
        candidate_y_mem = clCreateBuffer(m_gContext, CL_MEM_READ_ONLY, buffer_size_1D, NULL, NULL);
        if (LabelImageType::ImageDimension == 3) {
            candidate_z_mem = clCreateBuffer(m_gContext, CL_MEM_READ_ONLY, buffer_size_1D, NULL, NULL);
        }
        to_labels_mem = clCreateBuffer(m_gContext, CL_MEM_READ_ONLY, buffer_size_1D, NULL, NULL);

        // output buffer size
        curvature_flow_mem = clCreateBuffer(m_gContext, CL_MEM_READ_WRITE, float_buffer_size_1D, NULL, &err);
        testing_mem = clCreateBuffer(m_gContext, CL_MEM_READ_WRITE, sizeof (unsigned int) * aGPUInGlobalSize1D, NULL, &err);

        itk_handle_clEnqueueCreateBuffer(err);
        assert(err == CL_SUCCESS);

        err |= clEnqueueWriteBuffer(m_gCmd_queue, candidate_x_mem, CL_TRUE, 0, 
                buffer_size_1D, (void*) aGPUInCandidate[0], 0, NULL, NULL);
        itk_handle_clEnqueueWriteBuffer(err);
        assert(err == CL_SUCCESS);

        err |= clEnqueueWriteBuffer(m_gCmd_queue, candidate_y_mem, CL_TRUE, 0, 
                buffer_size_1D, (void*) aGPUInCandidate[1], 0, NULL, NULL);
        itk_handle_clEnqueueWriteBuffer(err);

        assert(err == CL_SUCCESS);
        if (LabelImageType::ImageDimension == 3) {
            err |= clEnqueueWriteBuffer(m_gCmd_queue, candidate_z_mem, CL_TRUE, 0, 
                    buffer_size_1D, (void*) aGPUInCandidate[2], 0, NULL, NULL);
            itk_handle_clEnqueueWriteBuffer(err);
        }

        err |= clEnqueueWriteBuffer(m_gCmd_queue, to_labels_mem, CL_TRUE, 0, 
                buffer_size_1D, (void*) aGPUInToLabel, 0, NULL, NULL);
        itk_handle_clEnqueueWriteBuffer(err);
        assert(err == CL_SUCCESS);

        size_t* origins = (size_t*) calloc(3, sizeof (size_t));

        size_t* region = (size_t*) calloc(3, sizeof (size_t));
        region[0] = aGPUInTotGlobSize[0];
        region[1] = aGPUInTotGlobSize[1];
        region[2] = aGPUInTotGlobSize[2]; // in 2D case this is 1.

        // the label image has been padded already, so the offset is 0. It could 
        // be potentially nice to not pad the image beforehands but to use an 
        // offset here. Although Eyad the master did not manage to do that 
        // (after a big effort).
        err |= clEnqueueWriteImage(m_gCmd_queue,
                labels_mem,
                CL_TRUE,
                origins,
                region,
                aGPUInTotGlobSize[0] * sizeof (int),
                (LabelImageType::ImageDimension == 2) ? 0 : slice * sizeof (int),
                (void*) aGPUInLabelImg,
                0, NULL, NULL);

        // error handling
        itk_handle_clEnqueueWriteImage(err);
        assert(err == CL_SUCCESS);

        // Finish everything in the queue before continuing
        err = clFinish(m_gCmd_queue);
        itk_handle_clFinish(err);
        assert(err == CL_SUCCESS);

        // What kernel to use?
        cl_kernel vGPUKernel = m_gKernels[0];

        // setup the arguments to our kernel
        unsigned int arg_num = 0;
        err = clSetKernelArg(vGPUKernel, arg_num++, sizeof (cl_mem), &labels_mem); // global memory (provide pointer)
        //        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (cl_mem), &masks_mem);  // global memeory
        //        err |= clSetKernelArg(vGPUKernel, arg_num++, mask_buffer_size, NULL);       // shared memory (no arg)

        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (cl_mem), &candidate_x_mem);
        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (cl_mem), &candidate_y_mem);
        if (LabelImageType::ImageDimension == 3) {
            err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (cl_mem), &candidate_z_mem);
        }
        err |= clSetKernelArg(vGPUKernel, arg_num++, local_buffer_size_1D, NULL);
        err |= clSetKernelArg(vGPUKernel, arg_num++, local_buffer_size_1D, NULL);
        if (LabelImageType::ImageDimension == 3) {
            err |= clSetKernelArg(vGPUKernel, arg_num++, local_buffer_size_1D, NULL);
        }

        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (cl_mem), &to_labels_mem);
        err |= clSetKernelArg(vGPUKernel, arg_num++, local_buffer_size_1D, NULL);

        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (float), &aGPUInVolume); // global memory
        //        err |= clSetKernelArg(vGPUKernel, arg_num++,  sizeof (float), &curvature_radius);
        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (float), &aGPUInEnergyContourLengthCoeff);

        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (float), &aGPUInMaskRadius[0]);
        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (float), &aGPUInMaskRadius[1]);
        if (LabelImageType::ImageDimension == 3) {
            err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (float), &aGPUInMaskRadius[2]);
        }

        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (unsigned int), &aGPUInOffset[0]);
        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (unsigned int), &aGPUInOffset[1]);
        if (LabelImageType::ImageDimension == 3) {
            err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (unsigned int), &aGPUInOffset[2]);
        }

        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (cl_mem), &curvature_flow_mem);
        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (cl_mem), &testing_mem);

        itk_handle_clSetKernelArgument(err);

        assert(err == CL_SUCCESS);


        size_t sizet_1D = aGPUInGlobalSize1D;
        size_t local_sizet_1D = aGPUInLocalSize1D;

        // Run the calculation by enqueuing it and forcing the command queue to complete the task

        err = clEnqueueNDRangeKernel(m_gCmd_queue,
                vGPUKernel, // kernel
                1, // dimension of execution model
                NULL,
                &sizet_1D, // padded problem size
                &local_sizet_1D, // work group size
                0, NULL, NULL);
        itk_handle_clEnqueueNDRangeKernel(err);

        // after this finshed, the calculations are for sure done
        clFinish(m_gCmd_queue);
        itk_handle_clFinish(err);
        assert(err == CL_SUCCESS);

        // Once finished read back the results from the answer
        // array into the results array
        err = clEnqueueReadBuffer(m_gCmd_queue, curvature_flow_mem, CL_TRUE, 0, float_buffer_size_1D,
                aGPUOutCurvatureFlow, 0, NULL, NULL);
        itk_handle_clEnqueueReadBuffer(err);
        assert(err == CL_SUCCESS);

        err = clEnqueueReadBuffer(m_gCmd_queue, testing_mem, CL_TRUE, 0, buffer_size_1D,
                aGPUOutTesting, 0, NULL, NULL);
        itk_handle_clEnqueueReadBuffer(err);
        assert(err == CL_SUCCESS);

        cl_int finish_err = clFinish(m_gCmd_queue);
        itk_handle_clFinish(err);
        assert(finish_err == CL_SUCCESS);

        err = clReleaseMemObject(labels_mem);
        err |= clReleaseMemObject(curvature_flow_mem);
        err |= clReleaseMemObject(testing_mem);
        err |= clReleaseMemObject(candidate_x_mem);
        err |= clReleaseMemObject(candidate_y_mem);
        if (LabelImageType::ImageDimension == 3) {
            err |= clReleaseMemObject(candidate_z_mem);
        }
        err |= clReleaseMemObject(to_labels_mem);
        itk_handle_clReleaseMemObject(err);
        assert(err == CL_SUCCESS);
        delete [] origins;
	delete [] region;

    }
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void RCGPUPiecewiseConstantSquareDistCurvatureRegEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::calculateExternalEnergyWithGPU(
    GPU_IntType* aGPUInLabelImg, float* aGPUInDataImg,
    unsigned int* aGPUInTotGlobSize,
    float aGPUInEnergyRegionCoeff, float aGPURegionMergingThreshold,
    float aGPUInBalloonForceCoeff, float aGPUInConstOutwardFlowCoeff,
    float* aGPUInMaskRadius, unsigned int* aGPUInOffset, float* aGPUOutEnergyDiff,
    unsigned int* aGPUOutMerging, unsigned int** aGPUInCandidate,
    unsigned int* aGPUInToLabel, unsigned int aGPUInGlobalSize1D, unsigned int aGPUInLocalSize,
    float* aGPUInMeans, float* aGPUInVariances, int aGPUInCounts) {

        cl_uint err;

        unsigned int slice = aGPUInTotGlobSize[0] * aGPUInTotGlobSize[1];

        //        size_t global_int_buffer_size_wg = sizeof (GPU_IntType) * aGPUInTotGlobSize[0] * aGPUInTotGlobSize[1];
        //        size_t global_float_buffer_size_wg = sizeof (GPU_IntType) * aGPUInTotGlobSize[0] * aGPUInTotGlobSize[1];
        size_t global_int_buffer_size_1D = sizeof (GPU_UIntType) * aGPUInGlobalSize1D;
        size_t global_uint_buffer_size_1D = sizeof (GPU_UIntType) * aGPUInGlobalSize1D;
        size_t global_float_buffer_size_1D = sizeof (float) * aGPUInGlobalSize1D;
        size_t local_int_buffer_size_1D = sizeof (GPU_UIntType) * aGPUInLocalSize; // gInLocalSize = workgroupsize
        size_t local_float_buffer_size_1D = sizeof (float) * aGPUInLocalSize;
        size_t label_count_float_buffer_size_1D = sizeof (float) * aGPUInCounts;
        //        size_t mask_buffer_size = sizeof (GPU_IntType) * (2 * aGPUInMaskRadius[0] + 1)*(2 * aGPUInMaskRadius[1] + 1);


        cl_mem labels_mem, data_mem, candidate_x_mem, candidate_y_mem, candidate_z_mem,
                to_labels_mem, mean_mem, var_mem,
                energyDiff_mem, merge_mem;



        cl_image_format format_int;
        format_int.image_channel_data_type = CL_SIGNED_INT32;
        format_int.image_channel_order = CL_R;

        cl_image_format format_float;
        format_float.image_channel_data_type = CL_FLOAT;
        format_float.image_channel_order = CL_R;

        cl_int imageerr;
        if (LabelImageType::ImageDimension == 2) {
            labels_mem = clCreateImage2D(
                    m_gContext, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                    &format_int,
                    aGPUInTotGlobSize[0], aGPUInTotGlobSize[1],
                    aGPUInTotGlobSize[0] * sizeof (int),
                    aGPUInLabelImg,
                    &imageerr);
            itk_handle_clEnqueueCreateImage(imageerr);

            data_mem = clCreateImage2D(
                    m_gContext,
                    CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                    &format_float,
                    aGPUInTotGlobSize[0], aGPUInTotGlobSize[1],
                    aGPUInTotGlobSize[0] * sizeof (float),
                    aGPUInDataImg,
                    &imageerr);
            itk_handle_clEnqueueCreateImage(imageerr);

        } else if (LabelImageType::ImageDimension == 3) {
            labels_mem = clCreateImage3D(m_gContext,
                    CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, &format_int,
                    aGPUInTotGlobSize[0], aGPUInTotGlobSize[1], aGPUInTotGlobSize[2],
                    aGPUInTotGlobSize[0] * sizeof (int),
                    slice * sizeof (int),
                    aGPUInLabelImg,
                    &imageerr);
            itk_handle_clEnqueueCreateImage(imageerr);

            data_mem = clCreateImage3D(m_gContext,
                    CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, &format_float,
                    aGPUInTotGlobSize[0], aGPUInTotGlobSize[1], aGPUInTotGlobSize[2],
                    aGPUInTotGlobSize[0] * sizeof (float),
                    slice * sizeof (float),
                    aGPUInDataImg,
                    &imageerr);
            itk_handle_clEnqueueCreateImage(imageerr);
        }
        // mask size
        //        masks_mem = clCreateBuffer(m_gContext, CL_MEM_READ_ONLY, mask_buffer_size, NULL, NULL);

        candidate_x_mem = clCreateBuffer(m_gContext, CL_MEM_READ_ONLY, global_int_buffer_size_1D, NULL, NULL);
        candidate_y_mem = clCreateBuffer(m_gContext, CL_MEM_READ_ONLY, global_int_buffer_size_1D, NULL, NULL);
        if (LabelImageType::ImageDimension == 3) {
            candidate_z_mem = clCreateBuffer(m_gContext, CL_MEM_READ_ONLY, global_int_buffer_size_1D, NULL, NULL);
        }
        to_labels_mem = clCreateBuffer(m_gContext, CL_MEM_READ_ONLY, global_int_buffer_size_1D, NULL, NULL);

        mean_mem = clCreateBuffer(m_gContext, CL_MEM_READ_ONLY, label_count_float_buffer_size_1D, NULL, NULL);
        var_mem = clCreateBuffer(m_gContext, CL_MEM_READ_ONLY, label_count_float_buffer_size_1D, NULL, NULL);

        // output buffer size
        energyDiff_mem = clCreateBuffer(m_gContext, CL_MEM_READ_WRITE, global_float_buffer_size_1D, NULL, NULL);
        merge_mem = clCreateBuffer(m_gContext, CL_MEM_READ_WRITE, global_int_buffer_size_1D, NULL, NULL);


        // creating buffers
        //        err = clEnqueueWriteBuffer(m_gCmd_queue, masks_mem, CL_TRUE, 0, mask_buffer_size, (void*) masks, 0, NULL, NULL);
        err = clEnqueueWriteBuffer(m_gCmd_queue, candidate_x_mem, CL_TRUE, 0, 
                global_int_buffer_size_1D, (void*) aGPUInCandidate[0], 0, NULL, NULL);
        err |= clEnqueueWriteBuffer(m_gCmd_queue, candidate_y_mem, CL_TRUE, 0, 
                global_int_buffer_size_1D, (void*) aGPUInCandidate[1], 0, NULL, NULL);
        if (LabelImageType::ImageDimension == 3) {
            err |= clEnqueueWriteBuffer(m_gCmd_queue, candidate_z_mem, CL_TRUE, 0, 
                    global_int_buffer_size_1D, (void*) aGPUInCandidate[2], 0, NULL, NULL);
        }
        err |= clEnqueueWriteBuffer(m_gCmd_queue, to_labels_mem, CL_TRUE, 0, 
                global_int_buffer_size_1D, (void*) aGPUInToLabel, 0, NULL, NULL);
        err |= clEnqueueWriteBuffer(m_gCmd_queue, mean_mem, CL_TRUE, 0, 
                label_count_float_buffer_size_1D, (void*) aGPUInMeans, 0, NULL, NULL);
        err |= clEnqueueWriteBuffer(m_gCmd_queue, var_mem, CL_TRUE, 0, 
                label_count_float_buffer_size_1D, (void*) aGPUInVariances, 0, NULL, NULL);
        itk_handle_clEnqueueWriteBuffer(err);

        size_t* origins = (size_t*) calloc(3, sizeof (size_t));
        //origins[0] = mask_radius[0];
        //origins[1] = mask_radius[1];
        size_t* region = (size_t*) calloc(3, sizeof (size_t));

        region[0] = aGPUInTotGlobSize[0];
        region[1] = aGPUInTotGlobSize[1];
        region[2] = aGPUInTotGlobSize[2]; // in 2D case this is 1.
        err = clEnqueueWriteImage(m_gCmd_queue,
                labels_mem,
                CL_TRUE,
                origins,
                region,
                aGPUInTotGlobSize[0] * sizeof (int),
                (LabelImageType::ImageDimension == 2) ? 0 : slice * sizeof (int),
                (void*) aGPUInLabelImg, 0, NULL, NULL);
        itk_handle_clEnqueueWriteImage(err);

        err = clEnqueueWriteImage(m_gCmd_queue,
                data_mem,
                CL_TRUE,
                origins,
                region,
                aGPUInTotGlobSize[0] * sizeof (float),
                (LabelImageType::ImageDimension == 2) ? 0 : slice * sizeof (float),
                (void*) aGPUInDataImg,
                0, NULL, NULL);
        itk_handle_clEnqueueWriteImage(err);
        assert(err == CL_SUCCESS);


        // Get all of the stuff written and allocated
        clFinish(m_gCmd_queue);
        itk_handle_clFinish(err);
        assert(err == CL_SUCCESS);

        // What kernel to use?
        cl_kernel vGPUKernel = m_gKernels[1];

        // Now setup the arguments to our kernel
        cl_uint arg_num = 0;
        err = clSetKernelArg(vGPUKernel, arg_num++, sizeof (cl_mem), &labels_mem);
        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (cl_mem), &data_mem);
        //        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (cl_mem), &masks_mem);
        //        err |= clSetKernelArg(vGPUKernel, arg_num++, mask_buffer_size, NULL);

        //        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (unsigned int), &gInSize1D);

        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (cl_mem), &candidate_x_mem);
        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (cl_mem), &candidate_y_mem);
        if (LabelImageType::ImageDimension == 3) {
            err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (cl_mem), &candidate_z_mem);
        }
        err |= clSetKernelArg(vGPUKernel, arg_num++, local_int_buffer_size_1D, NULL);
        err |= clSetKernelArg(vGPUKernel, arg_num++, local_int_buffer_size_1D, NULL);
        if (LabelImageType::ImageDimension == 3) {
            err |= clSetKernelArg(vGPUKernel, arg_num++, local_int_buffer_size_1D, NULL);
        }

        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (cl_mem), &to_labels_mem);
        err |= clSetKernelArg(vGPUKernel, arg_num++, local_int_buffer_size_1D, NULL);

        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (float), &aGPUInEnergyRegionCoeff);


        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (float), &aGPURegionMergingThreshold);

        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (float), &aGPUInBalloonForceCoeff);
        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (float), &aGPUInConstOutwardFlowCoeff);

        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (cl_mem), &mean_mem);
        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (cl_mem), &var_mem);

        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (float), &aGPUInMaskRadius[0]);
        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (float), &aGPUInMaskRadius[1]);
        if (LabelImageType::ImageDimension == 3) {
            err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (float), &aGPUInMaskRadius[2]);
        }

        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (unsigned int), &aGPUInOffset[0]);
        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (unsigned int), &aGPUInOffset[1]);
        if (LabelImageType::ImageDimension == 3) {
            err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (unsigned int), &aGPUInOffset[2]);
        }
        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (cl_mem), &energyDiff_mem);
        err |= clSetKernelArg(vGPUKernel, arg_num++, sizeof (cl_mem), &merge_mem);

        itk_handle_clSetKernelArgument(err);

        assert(err == CL_SUCCESS);


        size_t sizet_1D = aGPUInGlobalSize1D; // this is the (padded) number of particles
        size_t local_sizet_1D = aGPUInLocalSize; // this is the work group size

        // Run the calculation by enqueuing it and forcing the command queue to
        // complete the task
        err = clEnqueueNDRangeKernel(m_gCmd_queue, vGPUKernel, 1, NULL,
                &sizet_1D, &local_sizet_1D, 0, NULL, NULL);

        itk_handle_clEnqueueNDRangeKernel(err);

        clFinish(m_gCmd_queue);
        itk_handle_clFinish(err);
        assert(err == CL_SUCCESS);

        // Once finished read back the results from the answer
        // array into the results array
        err = clEnqueueReadBuffer(m_gCmd_queue, energyDiff_mem, CL_TRUE, 0, global_float_buffer_size_1D,
                aGPUOutEnergyDiff, 0, NULL, NULL);
        assert(err == CL_SUCCESS);

        err = clEnqueueReadBuffer(m_gCmd_queue, merge_mem, CL_TRUE, 0, global_int_buffer_size_1D,
                aGPUOutMerging, 0, NULL, NULL);
        assert(err == CL_SUCCESS);

        err = clFinish(m_gCmd_queue);
        itk_handle_clFinish(err);
        assert(err == CL_SUCCESS);

        err = clReleaseMemObject(labels_mem);
        err |= clReleaseMemObject(data_mem);
        err |= clReleaseMemObject(energyDiff_mem);
        err |= clReleaseMemObject(candidate_x_mem);
        err |= clReleaseMemObject(candidate_y_mem);
        if (LabelImageType::ImageDimension == 3) {
            err |= clReleaseMemObject(candidate_z_mem);
        }
        err |= clReleaseMemObject(to_labels_mem);
        err |= clReleaseMemObject(mean_mem);
        err |= clReleaseMemObject(var_mem);
        err |= clReleaseMemObject(merge_mem);

        delete [] origins;
	delete [] region;
    }

    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void RCGPUPiecewiseConstantSquareDistCurvatureRegEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::CleanUp() {

        delete [] m_gDataImage;
        delete [] m_gLabelImage;

        clReleaseCommandQueue(m_gCmd_queue);
        clReleaseContext(m_gContext);

    }
        


    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void RCGPUPiecewiseConstantSquareDistCurvatureRegEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::AddPoint(IndexType aInd, LabelAbsPixelType aLabel, DataPixelType aVal) {
        CountStatisticsIteratorType vIt = m_Count.find(aLabel);
        if(vIt != m_Count.end()) {
            vIt->second++;
            m_Sums[aLabel] += aVal;
            m_Sums_2[aLabel] += aVal * aVal;
        } else {
            m_Count[aLabel] = 1;
            m_Sums[aLabel] = aVal;
            m_Sums_2[aLabel] = aVal * aVal;
        }
    }     
    
    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void RCGPUPiecewiseConstantSquareDistCurvatureRegEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::RemovePoint(IndexType aInd /* = 0 */, LabelAbsPixelType aLabel, DataPixelType aVal) {
        assert(m_Count.find(aLabel) != m_Count.end());
        if(--m_Count[aLabel] == 0) {
            KillRegion(aLabel);
        } else {
            m_Sums[aLabel] -= aVal;
            m_Sums_2[aLabel] -= aVal * aVal;
        }
    }
    



    template<typename TLabelImage, typename TDataImage, typename TEnergyDifference>
    void RCGPUPiecewiseConstantSquareDistCurvatureRegEnergy<TLabelImage, TDataImage, TEnergyDifference>
    ::KillRegion(LabelAbsPixelType aLabel) {
        m_Count.erase(aLabel);
        m_Sums.erase(aLabel);
        m_Sums_2.erase(aLabel);
    }
        


} // end namespace itk
