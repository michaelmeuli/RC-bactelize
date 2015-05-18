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
float calculateKLMergingCriterion(float m1, float m2, float v1, float v2, int c1, int c2) {
  float m12 = (c1 * m1 + c2 * m2)/ (c1 + c2);
  float vs1 = v1 * (c1 - 1) + c1 * m1 * m1;
  float vs2 = v2 * (c2 - 1) + c2 * m2 * m2;

  float v12 = ( 1.0f / (c1 + c2 - 1.0f)) * (vs1 + vs2 - (c1 + c2) * m12 * m12);
  float kdl1 = (m1 - m12) * (m1 - m12) / (2.0f * v12) + 0.5f * (v1 / v12 - 1.0f - log(v1 / v12));
  float kdl2 = (m2 - m12) * (m2 - m12) / (2.0f * v12) + 0.5f * (v2 / v12 - 1.0f - log(v2 / v12));
	
  return kdl1 + kdl2;
}

__kernel void
energy(
       __read_only image3d_t label, // the label image
       __read_only image3d_t data,  // the data image

       __global uint *ram_candidate_x, // global particle coordinates 
       __global uint *ram_candidate_y,
       __global uint *ram_candidate_z,
	
       __local uint *candidate_x, // local particle coordinates
       __local uint *candidate_y,
       __local uint *candidate_z,

       __global uint *ram_to_label, // the candidate labels
       __local uint *to_label,      // local version of candidate labels

       float energy_region_coeff,      // data term coefficient
       float region_merging_threshold, // region merging ths 
       float balloon_force_coeff,      // balloon force coefficient
       float constant_outward_flow,    // outward flow coefficient 

       __global float *mean_label,     // an array of means (for each label)
       __global float *var_label,      // an array of variances (for each label)
	
       uint mask_radius_x,  // the mask radius (TODO:change to float)
       uint mask_radius_y,  // the offset is guaranteed to be >= radius
       uint mask_radius_z,

       uint offset0, // the offset from padding (ghost layer)
       uint offset1,
       uint offset2,

       __global float* energy_diff,  // energy output
       __global uint* merge        // local merging test (boolean)
       )
{	

  int gid = get_global_id(0);         // thread id (or particle id)
  int lid = get_local_id(0);          // id within the workgroup
  int group_size = get_local_size(0); // work-group size argument    

  candidate_x[lid] = ram_candidate_x[gid]; // coalesced copying to local mem
  candidate_y[lid] = ram_candidate_y[gid];
  candidate_z[lid] = ram_candidate_z[gid];
  to_label[lid] = ram_to_label[gid];
	
  sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;
	
  // check this one carefully since
  int4 coords = (int4)(candidate_x[lid] + offset0, candidate_y[lid] + offset1, candidate_z[lid] + offset2, 0);
  int4 values = read_imagei(label, sampler, coords);
  int from_label = abs(values.x);

  float pixel_data = (read_imagef(data, sampler, coords)).x;
	
  uint3 region_size = (uint3)(2 * mask_radius_x + 1, 2 * mask_radius_y + 1, 2 * mask_radius_z + 1);
	
  /*
    int number_times = (region_size.x * region_size.y * region_size.z + group_size - 1) / group_size;
    for (int i = 0;i < number_times;i++) {
    if ((lid + i * group_size) < region_size.x * region_size.y * region_size.z)
    masks[lid + i * group_size] = ram_masks[lid + i * group_size];
    barrier(CLK_LOCAL_MEM_FENCE);			
    }
  */


  merge[gid] = 0;
  int temp_from = 0;
  int temp_to = 0;
	
  int counter_from = -1;
  int counter_to = 0;
  float sum_from = -pixel_data;
  float sum_to = 0.0f;
  float double_sum_from = -pixel_data * pixel_data;
  float double_sum_to = 0.0f;
  for (int k = 0;k < region_size.z;k++) {
    for (int i = 0;i < region_size.y;i++) {
      for (int j = 0;j < region_size.x;j++) {
	int4 label_neighbor = read_imagei(label, sampler, 
					  (int4)(candidate_x[lid] + j + offset0 - mask_radius_x, 
						 candidate_y[lid] + i + offset1 - mask_radius_y, 
						 candidate_z[lid] + k + offset2 - mask_radius_z, 0));
	float data_neighbor = (read_imagef(data, sampler, 
					   (int4) (candidate_x[lid] + j + offset0 - mask_radius_x, 
						   candidate_y[lid] + i + offset1 - mask_radius_y, 
						   candidate_z[lid] + k + offset2 - mask_radius_z, 0))).x;

	// bank conflicts because of accessing arbitrary neighboring labels. look into that.      
	temp_from = from_label - abs(label_neighbor.x); 
	temp_to = to_label[lid] - abs(label_neighbor.x);

	if (temp_from != 0)
	  temp_from = 0;
	else
	  temp_from = 1;//masks[k * region_size.x * region_size.y + i * region_size.x + j];

	if (temp_to != 0)
	  temp_to = 0;
	else
	  temp_to = 1;//masks[k * region_size.x * region_size.y + i * region_size.x + j];

	float mask_val = 0;
	if((float)((k - mask_radius_z) * (k - mask_radius_z)) / (mask_radius_z * mask_radius_z) +
	   (float)((i - mask_radius_y) * (i - mask_radius_y)) / (mask_radius_y * mask_radius_y) +
	   (float)((j - mask_radius_x) * (j - mask_radius_x)) / (mask_radius_x * mask_radius_x) <= 1.0f)  {

	  mask_val = 1;
	}




	sum_from = sum_from + temp_from * data_neighbor * mask_val;
	sum_to = sum_to + temp_to * data_neighbor * mask_val;
	double_sum_from = double_sum_from + temp_from * data_neighbor * data_neighbor * mask_val;
	double_sum_to = double_sum_to + temp_to * data_neighbor * data_neighbor * mask_val;
	counter_from = counter_from + temp_from * mask_val; 
	counter_to = counter_to + temp_to * mask_val; 
      }
    }
  }

  float float_counter_from = (float)(counter_from);
  float float_counter_to = (float)(counter_to);

  float mean_from, var_from;
  float mean_to, var_to;

  if (counter_to == 0) {
    mean_to = mean_label[to_label[lid]];
    var_to = var_label[to_label[lid]];
  } else {
    mean_to = sum_to / float_counter_to;
    var_to = (double_sum_to - sum_to * sum_to / float_counter_to) / float_counter_to;
  }

  if (counter_from == 0) {
    mean_from = mean_label[from_label];
    var_from = var_label[from_label];
  } else {
    mean_from = sum_from / float_counter_from;
    var_from = (double_sum_from - sum_from * sum_from / float_counter_from) / float_counter_from;
  }

  energy_diff[gid] = energy_region_coeff * 
    ((pixel_data - mean_to) * (pixel_data - mean_to) - 
     (pixel_data - mean_from) * (pixel_data - mean_from));
  
  float criterion = calculateKLMergingCriterion(mean_from, 
						mean_to, 
						var_from, 
						var_to, 
						float_counter_from, 
						float_counter_to);

  merge[gid] = criterion < region_merging_threshold;
  if (from_label == 0) {
    energy_diff[gid] -= constant_outward_flow;
    if (balloon_force_coeff > 0)
      energy_diff[gid] -= (balloon_force_coeff * pixel_data);
    else
      energy_diff[gid] -= (balloon_force_coeff * (1 - pixel_data));
  } else {
    if (to_label[lid] == 0)
      energy_diff[gid] += constant_outward_flow;
  }

}

