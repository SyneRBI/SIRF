/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2019 University College London.

This is software developed for the Collaborative Computational
Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
(http://www.ccpsynerbi.ac.uk/).

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

#include "sirf/common/GeometricalInfo.h"
#include <iostream>
#include <sstream>
#include <string>

using namespace sirf;

template <int num_dimensions>
void
VoxelisedGeometricalInfo<num_dimensions>::
print_info() const
{
	std::cout << get_info();
}

template <int num_dimensions>
std::string
VoxelisedGeometricalInfo<num_dimensions>::
get_info() const
{
	std::stringstream info;
	info << "Offset: (";
	for (int i = 0; i < num_dimensions; i++) {
		info << _offset[i];
		if (i < num_dimensions - 1)
			info << ", ";
		else
			info << ")\n";
	}

	info << "Spacing: (";
	for (int i = 0; i < num_dimensions; i++) {
		info << _spacing[i];
		if (i < num_dimensions - 1)
			info << ", ";
		else
			info << ")\n";
	}

	info << "Size: (";
	for (int i = 0; i < num_dimensions; i++) {
		info << _size[i];
		if (i < num_dimensions - 1)
			info << ", ";
		else
			info << ")\n";
	}

	info << "Direction matrix: \n";
	for (int i = 0; i < num_dimensions; i++) {
		for (int j = 0; j < num_dimensions; j++) {
			info << _direction[i][j];
			if (j < num_dimensions - 1)
				info << ", ";
			else
				info << "\n";
		}
	}
	info << "\n";
	return info.str();
}

template <int num_dimensions>
VoxelisedGeometricalInfo<num_dimensions>::
VoxelisedGeometricalInfo(
	const Offset& offset, const Spacing& spacing,
	const Size& size, const DirectionMatrix& direction)
	:
	_offset(offset),
	_spacing(spacing),
	_size(size),
	_direction(direction)
{}

template <int num_dimensions>
const typename VoxelisedGeometricalInfo<num_dimensions>::Offset
VoxelisedGeometricalInfo<num_dimensions>::get_offset() const
{
	return _offset;
}

template <int num_dimensions>
const typename VoxelisedGeometricalInfo<num_dimensions>::Spacing
VoxelisedGeometricalInfo<num_dimensions>::get_spacing() const
{
	return _spacing;
}

template <int num_dimensions>
const typename VoxelisedGeometricalInfo<num_dimensions>::Size
VoxelisedGeometricalInfo<num_dimensions>::get_size() const
{
	return _size;
}

template <int num_dimensions>
const typename VoxelisedGeometricalInfo<num_dimensions>::DirectionMatrix
VoxelisedGeometricalInfo<num_dimensions>::get_direction() const
{
	return _direction;
}

template <int num_dimensions>
const typename VoxelisedGeometricalInfo<num_dimensions>::TransformMatrix
VoxelisedGeometricalInfo<num_dimensions>::
calculate_index_to_physical_point_matrix() const
{
	TransformMatrix index_to_physical_point_matrix;
	for (unsigned int dim = 0; dim<num_dimensions; dim++) {
		for (unsigned int axis = 0; axis<num_dimensions; axis++) {
			index_to_physical_point_matrix[dim][axis] =
				_direction[dim][axis] * _spacing[axis];
		}
		index_to_physical_point_matrix[dim][num_dimensions] = _offset[dim];
		index_to_physical_point_matrix[num_dimensions][dim] = 0;
	}
	index_to_physical_point_matrix[num_dimensions][num_dimensions] = 1;
	return index_to_physical_point_matrix;
}

namespace sirf {
template class VoxelisedGeometricalInfo<3>;
}
