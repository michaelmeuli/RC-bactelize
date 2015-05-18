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
 * Eyad Ebrahim
 *=========================================================================*/
//#pragma OPENCL EXTENSION cl_amd_printf : enable
// only shrinking. that means, the from_label is the abs(label). and to_label is '0'.
__kernel void
energy(
       __read_only image3d_t image,
       //	__global int* ram_masks,
       //	__local int* masks, // complies with y first , x second

       __global uint *ram_candidate_x,
       __global uint *ram_candidate_y,
       __global uint *ram_candidate_z,

       __local uint *candidate_x,
       __local uint *candidate_y,
       __local uint *candidate_z,

       __global uint *ram_to_label,
       __local uint *to_label,

       float volume,
       float curvature_mask_radius,
       float energy_contour_length_coeff,
	
       /** the max of both mask sizes */
       uint msize0, // y first, x second convention
       uint msize1,
       uint msize2,

       uint offset0,
       uint offset1,
       uint offset2,

       __global float* curvature_flow, // size of this one is different from the other kernels. it's the same size as the candidate lists.
       __global uint* testing)
{	

  int gid = get_global_id(0);

  int lid = get_local_id(0);
  int group_size = get_local_size(0);

  candidate_x[lid] = ram_candidate_x[gid];
  candidate_y[lid] = ram_candidate_y[gid];
  candidate_z[lid] = ram_candidate_z[gid];
  to_label[lid] = ram_to_label[gid];

  // maybe a barrier here. but I don't see the reason for that.
	
  sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;
	
  // check this one carefully since
  int4 coords = (int4)(candidate_x[lid] + offset0, candidate_y[lid] + offset1, candidate_z[lid] + offset2, 0);
  int4 values = read_imagei(image, sampler, coords);
  int from_label = abs(values.x); 

  uint3 region_size = (uint3)(2 * msize0 + 1, 2 * msize1 + 1, 2 * msize2 + 1);
	
  //	int number_times = (region_size.x * region_size.y * region_size.z + group_size - 1) / group_size;
  //	for (int i = 0;i < number_times;i++) {
  //		if ((lid + i * group_size) < region_size.x * region_size.y * region_size.z)
  //			masks[lid + i * group_size] = ram_masks[lid + i * group_size];

  //	}

  barrier(CLK_LOCAL_MEM_FENCE);
	
  testing[gid] = 0;

  int temp_from = 0;
  int temp_to = 0;
  int counter_from = 0;
  int counter_to = 0;

  for (int k = 0;k < region_size.z;k++) {
    for (int i = 0;i < region_size.y;i++) {
      for (int j = 0;j < region_size.x;j++) {

	// the accurate thing is the following
	int4 neighbor = read_imagei(image, sampler, 
				    (int4)(candidate_x[lid] + j + offset0 - msize0, 
					   candidate_y[lid] + i + offset1 - msize1, 
					   candidate_z[lid] + k + offset2 - msize2, 0));


	// next line: bank conflicts because of accessing arbitrary neighboring 
	// labels. look into that.
	temp_from = from_label - abs(neighbor.x); 
	temp_to = to_label[lid] - abs(neighbor.x);
	if (temp_from != 0)
	  temp_from = 0;
	else
	  temp_from = 1;

	if (temp_to != 0)
	  temp_to = 0;
	else
	  temp_to = 1;
	// counter_from = counter_from + temp_from * masks[k * region_size.x * region_size.y + i * region_size.x + j];
	// counter_to = counter_to + temp_to * masks[k * region_size.x * region_size.y + i * region_size.x + j];

	// float mask_val = (float)((k + msize2) * (k + msize2) // z coord
	// 			 + (i + msize0) * (i + msize0) // y coord 
	// 			 + (j + msize1) * (j + msize1) // z coord
	// 			 < curvature_mask_radius * curvature_mask_radius);

	float mask_val = 0;
	if((k + msize2) * (k + msize2) // z coord
	   + (i + msize0) * (i + msize0) // y coord 
	   + (j + msize1) * (j + msize1) // z coord
	   < curvature_mask_radius * curvature_mask_radius) 
	  mask_val = 1;
				  
	counter_from = counter_from + temp_from * mask_val;
	counter_to = counter_to + temp_to * mask_val;
      }
    }
  }
	
  float float_counter_from = (float)(counter_from);
  float float_counter_to = (float)(counter_to);

  if (from_label == 0)
    curvature_flow[gid] = -energy_contour_length_coeff * 16.0f / 
      (3.0f * curvature_mask_radius) * (float_counter_to / volume - 0.5f);
  else
    if (to_label[lid] == 0) {
      curvature_flow[gid] = energy_contour_length_coeff * 16.0f / 
	(3.0f * curvature_mask_radius) * (float_counter_from / volume - 0.5f);
    } else {
      curvature_flow[gid] = -energy_contour_length_coeff * 16.0f / 
	(3.0f * curvature_mask_radius) * (float_counter_to / volume - 0.5f);
      curvature_flow[gid] += energy_contour_length_coeff * 16.0f / 
	(3.0f * curvature_mask_radius) * (float_counter_from / volume - 0.5f);
    }
		
}
