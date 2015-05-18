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
__kernel void
energy(
       __read_only image3d_t image, // the label image
       //	__global int* ram_masks, // the mask of the effective neighborhood
       //	__local int* masks, // complies with y first , x second

       __global uint *ram_candidate_x,// global particle coordinates 
       __global uint *ram_candidate_y,
       __global uint *ram_candidate_z,

       __local uint *candidate_x, // local particle coordinates
       __local uint *candidate_y,
       __local uint *candidate_z,

       __global uint *ram_to_label, // the candidate labels
       __local uint *to_label, // local version of candidate labels

       float volume, // volume of the effective mask
       float energy_contour_length_coeff,
	
       uint mask_radius_x, // the mask radius (TODO:change to float)
       uint mask_radius_y,
       uint mask_radius_z,

       uint offset0, // the offset from padding (ghost layer)
       uint offset1, // the offset is guaranteed to be >= radius
       uint offset2,

       __global float* curvature_flow,  // energy output
       __global uint* testing // in and output for debugging
       )
{	

  int gid = get_global_id(0); // thread id (or particle id)

  int lid = get_local_id(0); // id within the workgroup
  int group_size = get_local_size(0); // work-group size argument

  candidate_x[lid] = ram_candidate_x[gid]; // coalesced copying to local mem
  candidate_y[lid] = ram_candidate_y[gid];
  candidate_z[lid] = ram_candidate_z[gid];
  to_label[lid] = ram_to_label[gid];
 
  // A sampler defines how to access images. Options:
  //  - no normalization: we don't want to have the indices between 0 and 1
  //  - CLK_ADDRESS_NONE: boundary management (we do it ourselves)
  //  - CLK_FILTER_NEAREST: no interpolation
  sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;
	
  // Read the label of the current particle from the label image
  int4 coords = (int4)(candidate_x[lid] + offset0, candidate_y[lid] + offset1, candidate_z[lid] + offset2, 0);
  int4 values = read_imagei(image, sampler, coords); // returns (x,y,z,w); y,z,w are empty here
  int from_label = abs(values.x); 

  uint3 region_size = (uint3)(2 * mask_radius_x + 1, 2 * mask_radius_y + 1, 2 * mask_radius_z + 1);

  // Copy the mask to local memory. This was replaced by compuatation since 
  // the local memory gets filled quite quickly for larger masks.	
  //	int number_times = (region_size.x * region_size.y * region_size.z + group_size - 1) / group_size;
  //	for (int i = 0;i < number_times;i++) {
  //		if ((lid + i * group_size) < region_size.x * region_size.y * region_size.z)
  //			masks[lid + i * group_size] = ram_masks[lid + i * group_size];
  //	}

  testing[gid] = 0; // to remove

  int temp_from = 0;
  int temp_to = 0;
  int counter_from = 0;
  int counter_to = 0;

  for (int k = 0;k < region_size.z;k++) {
    for (int i = 0;i < region_size.y;i++) {
      for (int j = 0;j < region_size.x;j++) {

	int4 neighbor = read_imagei(image, sampler, 
				    (int4)(candidate_x[lid] + j + offset0 - mask_radius_x, 
					   candidate_y[lid] + i + offset1 - mask_radius_y, 
					   candidate_z[lid] + k + offset2 - mask_radius_z, 0));


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
	// counter_from = counter_from + temp_from * 
	//   masks[k * region_size.x * region_size.y + i * region_size.x + j];
	// counter_to = counter_to + temp_to * 
	//   masks[k * region_size.x * region_size.y + i * region_size.x + j];

	// float mask_val = (float)((k + mask_radius_z) * (k + mask_radius_z) // z coord
	// 			 + (i + mask_radius_x) * (i + mask_radius_x) // y coord 
	// 			 + (j + mask_radius_y) * (j + mask_radius_y) // z coord
	// 			 < curvature_mask_radius * curvature_mask_radius);

	float mask_val = 0;
	if((float)((k - mask_radius_z) * (k - mask_radius_z)) / (mask_radius_z * mask_radius_z) +
	   (float)((i - mask_radius_y) * (i - mask_radius_y)) / (mask_radius_y * mask_radius_y) +
	   (float)((j - mask_radius_x) * (j - mask_radius_x)) / (mask_radius_x * mask_radius_x) <= 1.0f)  {

	  mask_val = 1;
	}
				  
	counter_from = counter_from + temp_from * mask_val;
	counter_to = counter_to + temp_to * mask_val;
      }
    }
  }
	
  // explicit casting to float
  float float_counter_from = (float)(counter_from);
  float float_counter_to = (float)(counter_to);

  // Calculate the energy values (see paper) and copy them to the global memory
  // As a radius (for normalization) we use the radius in x direction.
  if (from_label == 0)
    curvature_flow[gid] = -energy_contour_length_coeff * 16.0f / 
      (3.0f *  mask_radius_x) * (float_counter_to / volume - 0.5f);
  else
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
