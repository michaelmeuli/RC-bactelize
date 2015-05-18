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
#ifndef DIMENSION
#define DIMENSION 2
#endif
	
#if !(DIMENSION == 2 || DIMENSION == 3)
#error "OpenCL Kernel does only support DIMESNION==2 or DIMENSION==3."
#endif 	

__kernel void
energy(
#if DIMENSION == 2
       __read_only image2d_t image, // the label image
#elif DIMENSION == 3
       __read_only image3d_t image, // the label image
#endif
       //	__global int* ram_masks,
       //	__local int* masks, // complies with y first , x second

       __global uint *ram_candidate_x, // global particle coordinates 
       __global uint *ram_candidate_y,
#if DIMENSION==3
       __global uint *ram_candidate_z,
#endif

       __local uint *candidate_x,  // local particle coordinates
       __local uint *candidate_y,
#if DIMENSION==3
       __global uint *candidate_z,
#endif


       __global uint *ram_to_label,  // the candidate labels
       __local uint *to_label, // local version of candidate labels

       float volume,  // the volume of the effective mask (to replace)
       float energy_contour_length_coeff,
	
       uint mask_radius_x, // the mask radius (TODO:change to float)
       uint mask_radius_y,
#if DIMENSION==3
       uint mask_radius_z,
#endif

       uint offset0,  // the offset from padding (ghost layer)
       uint offset1,
#if DIMENSION==3
       uint offset2,
#endif

       __global float* curvature_flow, // the energy output
       __global uint* testing
       )         
{	

  int gid = get_global_id(0);         // thread id (or particle id)
  int lid = get_local_id(0);          // id within the workgroup
  int group_size = get_local_size(0); // work-group size

  candidate_x[lid] = ram_candidate_x[gid];
  candidate_y[lid] = ram_candidate_y[gid];
#if DIMENSION==3
  candidate_z[lid] = ram_candidate_z[gid];
#endif

  to_label[lid] = ram_to_label[gid];

  sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;

#if DIMENSION == 2
  int2 coords = (int2)(candidate_x[lid] + offset0, candidate_y[lid] + offset1);
  uint2 region_size = (uint2)(2 * mask_radius_x + 1, 2 * mask_radius_y + 1);
#elif DIMENSION == 3
  int4 coords = (int4)(candidate_x[lid] + offset0, candidate_y[lid] + offset1, candidate_z[lid] + offset2, 0); 
  uint3 region_size = (uint3)(2 * mask_radius_x + 1, 2 * mask_radius_y + 1, 2 * mask_radius_z + 1);
#endif

  int4 values = read_imagei(image, sampler, coords);
  int from_label = abs(values.x); 

  testing[gid] = 0;

  int temp_from = 0;
  int temp_to = 0;
  int counter_from = 0;
  int counter_to = 0;

#if DIMENSION==3
  for (int k = 0;k < region_size.z;k++) {
#endif
    for (int i = 0;i < region_size.y;i++) {
      for (int j = 0;j < region_size.x;j++) {

#if DIMENSION==2
	int2 neigh_coords = (int2)(candidate_x[lid] + j + offset0 - mask_radius_x,
				   candidate_y[lid] + i + offset1 - mask_radius_y);
#elif DIMENSION == 3
	int4 neigh_coords = (int4)(candidate_x[lid] + j + offset0 - mask_radius_x,
				   candidate_y[lid] + i + offset1 - mask_radius_y,
				   candidate_z[lid] + k + offset2 - mask_radius_z, 0);
#endif
	int4 neighbor = read_imagei(image, sampler, neigh_coords);
	
	// next line: potential bank conflicts because of accessing arbitrary neighboring labels. 
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

	// squared 'radius' of the ellipsoid
	float r_2 = 0.0f;
	r_2 += (float)((i - mask_radius_y) * (i - mask_radius_y)) / (mask_radius_y * mask_radius_y);
	r_2 += (float)((j - mask_radius_x) * (j - mask_radius_x)) / (mask_radius_x * mask_radius_x);
#if DIMENSION == 3
	r_2 += (float)((k - mask_radius_z) * (k - mask_radius_z)) / (mask_radius_z * mask_radius_z);
#endif
	float mask_val = 0;
	if( r_2 <= 1.0f)  { 
	  mask_val = 1;
	}
	
	counter_from = counter_from + temp_from * mask_val;
	counter_to = counter_to + temp_to * mask_val;
      }
    }
#if DIMENSION == 3
  }
#endif
  
  float float_counter_to = (float)(counter_to);
  float float_counter_from = (float)(counter_from);


  // Calculate the energy values (see paper) and copy them to the global memory
  // As a radius (for normalization) we use the radius in x direction.
#if DIMENSION == 2

  if (from_label == 0) {
    curvature_flow[gid] = -energy_contour_length_coeff * 3.0f * 3.141592f / 
      mask_radius_x * (float_counter_to / volume - 0.5f);
  } else {
    if (to_label[lid] == 0) {
      curvature_flow[gid] = energy_contour_length_coeff * 3.0f * 3.141592f / 
	mask_radius_x * (float_counter_from / volume - 0.5f);
    } else {
      curvature_flow[gid] = -energy_contour_length_coeff * 3.0f * 3.141592f / 
	mask_radius_x * (float_counter_to / volume - 0.5f);
      curvature_flow[gid] += energy_contour_length_coeff * 3.0f * 3.141592f / 
	mask_radius_x * (float_counter_from / volume - 0.5f);
    }
  }
#elif DIMENSION == 3

  if (from_label == 0) {
    curvature_flow[gid] = -energy_contour_length_coeff * 16.0f / 
      (3.0f *  mask_radius_x) * (float_counter_to / volume - 0.5f);
  } else {
    if (to_label[lid] == 0) {
      curvature_flow[gid] = energy_contour_length_coeff * 16.0f / 
	(3.0f * mask_radius_x) * (float_counter_from / volume - 0.5f);
    } else {
      curvature_flow[gid] = -energy_contour_length_coeff * 16.0f / 
	(3.0f * mask_radius_x ) * (float_counter_to / volume - 0.5f);
      curvature_flow[gid] += energy_contour_length_coeff * 16.0f / 
	(3.0f * mask_radius_x) * (float_counter_from / volume - 0.5f);
    }
  }
#endif
}
