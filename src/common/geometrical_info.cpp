/*
  CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
  Copyright 2018 Commonwealth Scientific and Industrial Research Organisation's
  Australian eHealth Research Organisation

  This is software developed for the Collaborative Computational
  Project in Positron Emission Tomography and Magnetic Resonance imaging
  (http://www.ccppetmr.ac.uk/).

  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
*/

#include <string>
#include "sirf/common/geometrical_info.h"

//template <int num_dimensions>
//VoxelisedGeometricalInfo<num_dimensions>::
//VoxelisedGeometricalInfo(
//	const Offset& offset, const Spacing& spacing,
//	const Size& size, const DirectionMatrix& direction)
//	:
//	offset (offset),
//	spacing (spacing),
//	size (size),
//	direction (direction)
//{}
//
//template <int num_dimensions>
//const typename VoxelisedGeometricalInfo<num_dimensions>::Offset
//VoxelisedGeometricalInfo<num_dimensions>::get_offset()
//{
//	return offset;
//}
//
//template <int num_dimensions>
//const typename VoxelisedGeometricalInfo<num_dimensions>::Spacing
//VoxelisedGeometricalInfo<num_dimensions>::get_spacing()
//{
//	return spacing;
//}
//
//template <int num_dimensions>
//const typename VoxelisedGeometricalInfo<num_dimensions>::Size
//VoxelisedGeometricalInfo<num_dimensions>::get_size()
//{
//	return size;
//}
//
//template <int num_dimensions>
//const typename VoxelisedGeometricalInfo<num_dimensions>::DirectionMatrix
//VoxelisedGeometricalInfo<num_dimensions>::get_direction()
//{
//	return direction;
//}
//
//template <int num_dimensions>
//const typename VoxelisedGeometricalInfo<num_dimensions>::TransformMatrix
//VoxelisedGeometricalInfo<num_dimensions>::
//calculate_index_to_physical_point_matrix()
//{
//	TransformMatrix index_to_physical_point_matrix;
//	for (unsigned int j=0; j<num_dimensions; j++) {
//		for (unsigned int i=0; i<num_dimensions; i++) {
//			// set cosines
//			index_to_physical_point_matrix[j][i] =
//				direction[j][i] * spacing[j];
//		}
//		// set translations
//		index_to_physical_point_matrix[j][num_dimensions] = offset[j];
//	}
//	index_to_physical_point_matrix[num_dimensions][num_dimensions] = 1;
//	return index_to_physical_point_matrix;
//}
//
