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
 * Eyad Ebrahim, Janick Cardinale
 *=========================================================================*/
#ifndef itkOpenCLSupporter_h
#define itkOpenCLSupporter_h

#include <stdio.h>
#include <assert.h>
#include <sys/sysctl.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>


#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

namespace itk {
    char * itk_clLoadProgramSource(const char *filename);

    int itk_init_clPlatform(cl_platform_id* platform);

    int itk_init_clDevice(cl_platform_id platform, cl_device_id& device);

    int itk_init_clContextAndQueue(cl_platform_id platform, cl_device_id device, cl_context& context,
                cl_command_queue& cmd_queue);

    int itk_create_clKernel(cl_context context, cl_kernel* kernel, cl_program* program, cl_device_id device,
            const char* filename, const char * kernelname);

    int itk_create_clKernel(cl_context context, cl_kernel* kernel, cl_program* program, cl_device_id device,
            const char* filename, const char * kernelname, unsigned int dimension, int index);

    int itk_handle_clEnqueueNDRangeKernel(int err);
    
    int itk_handle_clEnqueueWriteBuffer(int err);
    int itk_handle_clEnqueueReadBuffer(int err); 
    
    int itk_handle_clEnqueueCreateBuffer(int err);

    int itk_handle_clEnqueueWriteImage(int err);

    int itk_handle_clEnqueueCreateImage(int err);

    int itk_handle_clSetKernelArgument(int err);

    int itk_handle_clFinish(int err);

    int itk_handle_clReleaseMemObject(int err);
    
}

#include "itkOpenCLSupporter.hxx"


#endif // itkOpenCLSupporter_h