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
#ifndef openclInitialization_hxx
#define openclInitialization_hxx

#include "itkOpenCLSupporter.h"


namespace itk{

char * itk_clLoadProgramSource(const char *filename)
{
	struct stat statbuf;
	FILE *fh;
	char *source;

	fh = fopen(filename, "r");
	if (fh == 0)
		return 0;

	stat(filename, &statbuf);
	source = (char *) malloc(statbuf.st_size + 1);
	size_t ignore_value = fread(source, statbuf.st_size, 1, fh);
	source[statbuf.st_size] = '\0';

	return source;
}




int itk_init_clPlatform(cl_platform_id* platform) {

    cl_uint numPlatforms = 1;
    (*platform) = NULL;
    cl_int status = clGetPlatformIDs(1, NULL, &numPlatforms);
    
    if (0 < numPlatforms)
    {
        cl_platform_id* platforms = (cl_platform_id*)malloc(sizeof(cl_platform_id) * numPlatforms);
        status = clGetPlatformIDs(numPlatforms, platforms, &numPlatforms);

        printf("Number of platforms found: %d\n", numPlatforms);
        for (unsigned i = 0; i < numPlatforms; ++i)
        {
            char pbuf[100];
            status = clGetPlatformInfo(platforms[i],
                                       CL_PLATFORM_VENDOR,
                                       sizeof(pbuf),
                                       pbuf,
                                       NULL);

            assert(status==CL_SUCCESS);

            
            *platform = platforms[i];
            printf("%s\n", pbuf);
            if (!strcmp(pbuf, "Advanced Micro Devices, Inc."))
            {
                break;
            }
        }
        free(platforms);
    }

    if( NULL == *platform)
    {
        printf("NULL platform found so Exiting Application.");
	return -1;
    }

    return 0;
}

/**
 * Device is only GPU, no cpu. In case of no gpu exists EXCEPTION.
 */

int itk_init_clDevice(cl_platform_id platform, cl_device_id* device_ptr) {
	cl_int err = 0;
        size_t returned_size;
	// Find the GPU CL device, this is what we really want
	// If there is no GPU device is CL capable, fall back to CPU
	err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, device_ptr, NULL);

        printf("error in finding a gpu: %d and the other error of no gpu found: %d\n", err , CL_DEVICE_NOT_FOUND);


	printf("CL_INVALID_PLATFORM: %d\n", CL_INVALID_PLATFORM);
	printf("CL_INVALID_DEVICE_TYPE: %d\n", CL_INVALID_DEVICE_TYPE);
	printf("CL_INVALID_VALUE: %d\n", CL_INVALID_VALUE);
	printf("CL_DEVICE_NOT_FOUND: %d\n", CL_DEVICE_NOT_FOUND);

	if (err != CL_SUCCESS)
		return -1;
	assert(*device_ptr);

	// Get some information about the returned device
	cl_char vendor_name[1024] = {0};
	cl_char device_name[1024] = {0};
	err = clGetDeviceInfo(*device_ptr, CL_DEVICE_VENDOR, sizeof(vendor_name),
						  vendor_name, &returned_size);
	err |= clGetDeviceInfo(*device_ptr, CL_DEVICE_NAME, sizeof(device_name),
						  device_name, &returned_size);
	assert(err == CL_SUCCESS);
	printf("Connecting to %s %s...\n", vendor_name, device_name);

	return 0;
}



int itk_init_clContextAndQueue(cl_platform_id platform, cl_device_id device, cl_context* context,
            cl_command_queue* cmd_queue)
{

        /*
         * If we could find our platform, use it.
         * Otherwise use just available platform.
         */

        cl_context_properties cps[3] =
        {
            CL_CONTEXT_PLATFORM,
            (cl_context_properties)platform,
            0
        };

	cl_int err = 0;
	// Now create a context to perform our calculation with the
	// specified device
	*context = clCreateContext(cps, 1, &device, NULL, NULL, &err);
	assert(err == CL_SUCCESS);

	// And also a command queue for the context
	*cmd_queue = clCreateCommandQueue(*context, device, 0, NULL);
	return 0;
}


int itk_create_clKernel(cl_context context, cl_kernel* kernel, cl_program* program,
            cl_device_id device, const char* filename, const char * kernelname)
{

    size_t log_len;
    char* build_log;
    cl_int err = 0;
    // Load the program source from disk
    // The kernel/program is the project directory and in Xcode the executable
    // is set to launch from that directory hence we use a relative path
    char *program_source = itk_clLoadProgramSource(filename);
    program[0] = clCreateProgramWithSource(context, 1, (const char**) &program_source,
            NULL, &err);
    assert(err == CL_SUCCESS);

    err = clBuildProgram(program[0], 0, NULL, NULL, NULL, NULL);
    printf("Error in building program: %d\n", err);

    err = clGetProgramBuildInfo(program[0], device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_len);

    build_log = (char*) malloc(sizeof (char) *(log_len + 1));
    err = clGetProgramBuildInfo(program[0], device, CL_PROGRAM_BUILD_LOG, log_len, build_log, NULL);

    printf("Error in building kernal: %s\n", build_log);
    assert(err == CL_SUCCESS);

    // Now create the kernel "objects" that we want to use in the example file
    kernel[0] = clCreateKernel(program[0], kernelname, &err);
    free(program_source);
    return 0;
}

int itk_create_clKernel(cl_context context, cl_kernel* kernel, cl_program* program,
            cl_device_id device, const char* filename, const char * kernelname,
            unsigned int dimension, int index)
{

    size_t log_len;
    char* build_log;
    cl_int err = 0;
    // Load the program source from disk
    // The kernel/program is the project directory and in Xcode the executable
    // is set to launch from that directory hence we use a relative path
    char *program_source = itk_clLoadProgramSource(filename);
    program[index] = clCreateProgramWithSource(context, 1, (const char**) &program_source,
            NULL, &err);
    assert(err == CL_SUCCESS);

    std::stringstream vSS;
    vSS << "-D DIMENSION=" << dimension;

    err = clBuildProgram(program[index], 0, NULL,
            const_cast<char *> (vSS.str().c_str()), NULL, NULL);
    printf("Error in building program: %d\n", err);

    err = clGetProgramBuildInfo(program[index], device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_len);

    build_log = (char*) malloc(sizeof (char) *(log_len + 1));
    err = clGetProgramBuildInfo(program[index], device, CL_PROGRAM_BUILD_LOG, log_len, build_log, NULL);

    printf("Error in building kernel: %s\n", build_log);
    assert(err == CL_SUCCESS);

    // create the kernel "objects" that we want to use
    kernel[index] = clCreateKernel(program[index], kernelname, &err);
    free(program_source);
    return 0;
}

int itk_handle_clEnqueueNDRangeKernel(int err) {
    if (err == CL_SUCCESS)
        return 0;
    printf("Error received from Enqueue ND Range Kernel: %d. Check the following table to compare:\n", err);

    printf("CL_INVALID_PROGRAM_EXECUTABLE: %d\n", CL_INVALID_PROGRAM_EXECUTABLE);
    printf("CL_INVALID_COMMAND_QUEUE: %d\n", CL_INVALID_COMMAND_QUEUE);
    printf("CL_INVALID_KERNEL: %d\n", CL_INVALID_KERNEL);
    printf("CL_INVALID_CONTEXT: %d\n", CL_INVALID_CONTEXT);
    printf("CL_INVALID_KERNEL_ARGS: %d\n", CL_INVALID_KERNEL_ARGS);
    printf("CL_INVALID_WORK_DIMENSION: %d\n", CL_INVALID_WORK_DIMENSION);
    printf("CL_INVALID_WORK_GROUP_SIZE: %d\n", CL_INVALID_WORK_GROUP_SIZE);
    printf("CL_INVALID_WORK_ITEM_SIZE: %d\n", CL_INVALID_WORK_ITEM_SIZE);
    printf("CL_INVALID_GLOBAL_OFFSET: %d\n", CL_INVALID_GLOBAL_OFFSET);
    printf("CL_OUT_OF_RESOURCES: %d\n", CL_OUT_OF_RESOURCES);
    printf("CL_MEM_OBJECT_ALLOCATION_FAILURE: %d\n", CL_MEM_OBJECT_ALLOCATION_FAILURE);
    printf("CL_INVALID_EVENT_WAIT_LIST: %d\n", CL_INVALID_EVENT_WAIT_LIST);
    printf("CL_OUT_OF_HOST_MEMORY: %d\n", CL_OUT_OF_HOST_MEMORY);

//    assert(err == CL_SUCCESS);

    return 0;
}

int itk_handle_clEnqueueWriteBuffer(int err) {
    if (err == CL_SUCCESS)
        return 0;
    printf("Error received from Enqueue Write Buffer: %d. Check the following table to compare:\n", err);

    printf("CL_INVALID_COMMAND_QUEUE: %d\n", CL_INVALID_COMMAND_QUEUE);
    printf("CL_INVALID_CONTEXT: %d\n", CL_INVALID_CONTEXT);
    printf("CL_INVALID_MEM_OBJECT: %d\n", CL_INVALID_MEM_OBJECT);
    printf("CL_INVALID_VALUE: %d\n", CL_INVALID_VALUE);
    printf("CL_INVALID_EVENT_WAIT_LIST: %d\n", CL_INVALID_EVENT_WAIT_LIST);
    printf("CL_MEM_OBJECT_ALLOCATION_FAILURE: %d\n", CL_MEM_OBJECT_ALLOCATION_FAILURE);
    printf("CL_OUT_OF_HOST_MEMORY: %d\n", CL_OUT_OF_HOST_MEMORY);

//    assert(err == CL_SUCCESS);

    return 0;
}

int itk_handle_clEnqueueReadBuffer(int err) {
    if (err == CL_SUCCESS)
        return 0;
    printf("Error received from Enqueue Read Buffer: %d. Check the following table to compare:\n", err);

    printf("CL_INVALID_COMMAND_QUEUE: %d\n", CL_INVALID_COMMAND_QUEUE);
    printf("CL_INVALID_CONTEXT: %d\n", CL_INVALID_CONTEXT);
    printf("CL_INVALID_MEM_OBJECT: %d\n", CL_INVALID_MEM_OBJECT);
    printf("CL_INVALID_VALUE: %d\n", CL_INVALID_VALUE);
    printf("CL_INVALID_EVENT_WAIT_LIST: %d\n", CL_INVALID_EVENT_WAIT_LIST);
    printf("CL_MEM_OBJECT_ALLOCATION_FAILURE: %d\n", CL_MEM_OBJECT_ALLOCATION_FAILURE);
    printf("CL_OUT_OF_RESOURCES: %d\n", CL_OUT_OF_RESOURCES);
    printf("CL_OUT_OF_HOST_MEMORY: %d\n", CL_OUT_OF_HOST_MEMORY);

//    assert(err == CL_SUCCESS);

    return 0;
}

int itk_handle_clEnqueueCreateBuffer(int err) {
    if (err == CL_SUCCESS)
        return 0;
    printf("Error received from Enqueue Create Buffer: %d. Check the following table to compare:\n", err);

    printf("CL_INVALID_CONTEXT: %d\n", CL_INVALID_CONTEXT);
    printf("CL_INVALID_BUFFER_SIZE: %d\n", CL_INVALID_BUFFER_SIZE);
    printf("CL_INVALID_MEM_OBJECT: %d\n", CL_INVALID_MEM_OBJECT);
    printf("CL_INVALID_VALUE: %d\n", CL_INVALID_VALUE);
    printf("CL_INVALID_HOST_PTR: %d\n", CL_INVALID_HOST_PTR);
    printf("CL_MEM_OBJECT_ALLOCATION_FAILURE: %d\n", CL_MEM_OBJECT_ALLOCATION_FAILURE);
    printf("CL_OUT_OF_HOST_MEMORY: %d\n", CL_OUT_OF_HOST_MEMORY);

//    assert(err == CL_SUCCESS);

    return 0;
}


int itk_handle_clEnqueueWriteImage(int err) {
    if (err == CL_SUCCESS)
        return 0;
    printf("Error received from Enqueue Write Image: %d. Check the following table to compare:\n", err);

    printf("CL_INVALID_CONTEXT: %d\n", CL_INVALID_CONTEXT);
    printf("CL_INVALID_MEM_OBJECT: %d\n", CL_INVALID_MEM_OBJECT);
    printf("CL_INVALID_VALUE: %d\n", CL_INVALID_VALUE);
    printf("CL_MEM_OBJECT_ALLOCATION_FAILURE: %d\n", CL_MEM_OBJECT_ALLOCATION_FAILURE);
    printf("CL_OUT_OF_HOST_MEMORY: %d\n", CL_OUT_OF_HOST_MEMORY);
    printf("CL_INVALID_COMMAND_QUEUE: %d\n", CL_INVALID_COMMAND_QUEUE);
    printf("CL_INVALID_EVENT_WAIT_LIST: %d\n", CL_INVALID_EVENT_WAIT_LIST);


//    assert(err == CL_SUCCESS);

    return 0;
}

int itk_handle_clEnqueueReadImage(int err) {
  if (err == CL_SUCCESS) {
        return 0;
  }
    printf("Error received from Enqueue Read Image: %d. Check the following table to compare:\n", err);

    printf("CL_INVALID_COMMAND_QUEUE: %d\n", CL_INVALID_COMMAND_QUEUE);
    printf("CL_INVALID_CONTEXT : %d\n", CL_INVALID_CONTEXT);
    printf("CL_INVALID_MEM_OBJECT : %d\n", CL_INVALID_MEM_OBJECT);
    printf("CL_INVALID_VALUE : %d\n", CL_INVALID_VALUE);
    printf("CL_INVALID_EVENT_WAIT_LIST : %d\n", CL_INVALID_EVENT_WAIT_LIST);
    printf("CL_MEM_OBJECT_ALLOCATION_FAILURE : %d\n", CL_MEM_OBJECT_ALLOCATION_FAILURE);
    printf("CL_OUT_OF_HOST_MEMORY : %d\n", CL_OUT_OF_HOST_MEMORY);

    return 0;
}


int itk_handle_clEnqueueCreateImage(int err) {
  if (err == CL_SUCCESS) {
        return 0;
  }
    printf("Error received from Enqueue Create Image: %d. Check the following table to compare:\n", err);

    printf("CL_INVALID_CONTEXT: %d\n", CL_INVALID_CONTEXT);
    printf("CL_INVALID_VALUE: %d\n", CL_INVALID_VALUE);
    printf("CL_INVALID_IMAGE_FORMAT_DESCRIPTOR: %d\n", CL_INVALID_IMAGE_FORMAT_DESCRIPTOR);
    printf("CL_INVALID_IMAGE_SIZE: %d\n", CL_INVALID_IMAGE_SIZE);
    printf("CL_INVALID_HOST_PTR: %d\n", CL_INVALID_HOST_PTR);
    printf("CL_IMAGE_FORMAT_NOT_SUPPORTED: %d\n", CL_IMAGE_FORMAT_NOT_SUPPORTED);
    printf("CL_INVALID_OPERATION: %d\n", CL_INVALID_OPERATION);
    printf("CL_MEM_OBJECT_ALLOCATION_FAILURE: %d\n", CL_MEM_OBJECT_ALLOCATION_FAILURE);
    printf("CL_OUT_OF_HOST_MEMORY: %d\n", CL_OUT_OF_HOST_MEMORY);


//    assert(err == CL_SUCCESS);

    return 0;
}


int itk_handle_clSetKernelArgument(int err) {
    if (err == CL_SUCCESS)
        return 0;
    printf("Error received from Enqueue Set Kernel Argument: %d. Check the following table to compare:\n", err);

    printf("CL_INVALID_KERNEL: %d\n", CL_INVALID_KERNEL);
    printf("CL_INVALID_ARG_INDEX: %d\n", CL_INVALID_ARG_INDEX);
    printf("CL_INVALID_ARG_VALUE: %d\n", CL_INVALID_ARG_VALUE);
    printf("CL_INVALID_MEM_OBJECT: %d\n", CL_INVALID_MEM_OBJECT);
    printf("CL_INVALID_SAMPLER: %d\n", CL_INVALID_SAMPLER);
    printf("CL_INVALID_ARG_SIZE: %d\n", CL_INVALID_ARG_SIZE);

//    assert(err == CL_SUCCESS);
    
    return 0;
}

int itk_handle_clFinish(int err) {
    if (err == CL_SUCCESS)
        return 0;
    printf("Error received from clFinish: %d. Check the following table to compare:\n", err);

    printf("CL_INVALID_COMMAND_QUEUE: %d\n", CL_INVALID_COMMAND_QUEUE);
    printf("CL_OUT_OF_HOST_MEMORY: %d\n", CL_OUT_OF_HOST_MEMORY);
    printf("CL_OUT_OF_RESOURCES: %d\n", CL_OUT_OF_RESOURCES);
//    assert(err == CL_SUCCESS);

    return 0;
}

int itk_handle_clReleaseMemObject(int err){
    if (err == CL_SUCCESS)
        return 0;
    printf("Error received from clReleaseMemObject: %d. Check the following table to compare:\n", err);
    printf("CL_INVALID_MEM_OBJECT: %d\n", CL_INVALID_MEM_OBJECT);

    return 0;
}

} // namespace itk
#endif // openclInitialization_hxx
