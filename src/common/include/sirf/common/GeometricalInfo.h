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

#ifndef SIRF_GEOMETRICAL_INFO_TYPE
#define SIRF_GEOMETRICAL_INFO_TYPE

namespace sirf {

template <int num_dimensions, typename T>
class tVector {
private:
	typedef T _vectT[num_dimensions];

public:
	inline T& operator[](const int d)
		{ return vect[d]; }
	inline T operator[](const int d) const
		{ return vect[d]; }

private:
	_vectT vect;
};


template <int num_physical_dimensions, int num_index_dimensions>
class GeometricalInfo {
public:
	typedef tVector<num_physical_dimensions, float>     Coordinate;
	typedef tVector<num_index_dimensions, unsigned int> Index;
	// Eventually something here like
	// Coordinate transform_index_to_physical_point(Index)
	// Index transform_physical_point_to_index(Coordinate)
};


template <int num_dimensions>
class VoxelisedGeometricalInfo :
	public GeometricalInfo<num_dimensions, num_dimensions> {
private:
	typedef GeometricalInfo<num_dimensions, num_dimensions> BaseType;

public:
	// TODO: Why to I have to define these again?
	typedef typename BaseType::Coordinate Coordinate;
	typedef typename BaseType::Index Index;

	/*!
	   Offset is the coordinate of the center of the first voxel in physical
	   space.
	*/
	typedef Coordinate Offset;
	/*!
	   Spacing is the physical distance between voxels in each dimensions.
	*/
	typedef Coordinate Spacing;
	/*!
	   Size is the number of voxels in each dimension.
	*/
	typedef Index Size;
	/*!
	   Each vector in Direction tells the direction of the axis in LPS
	   physical space.
	*/
	typedef tVector<num_dimensions, Coordinate> DirectionMatrix;
	/*!
	   Each vector in Direction tells the direction of the axis in LPS
	   physical space.
	*/
	typedef tVector<num_dimensions+1, tVector<num_dimensions+1, float> >
		TransformMatrix;


	VoxelisedGeometricalInfo(
		const Offset& offset, const Spacing& spacing,
		const Size& size, const DirectionMatrix& direction);
	virtual ~VoxelisedGeometricalInfo() {};

    const Offset get_offset() const;
    const Spacing get_spacing() const;
    const Size get_size() const;
    const DirectionMatrix get_direction() const;

    const TransformMatrix calculate_index_to_physical_point_matrix() const;

private:
	Offset offset;
	Spacing spacing;
	Size size;
	DirectionMatrix direction;
};

typedef GeometricalInfo<3, 3> GeometricalInfo3D;
typedef VoxelisedGeometricalInfo<3> VoxelisedGeometricalInfo3D;
typedef VoxelisedGeometricalInfo<3>::TransformMatrix TransformMatrix3D;

template <int num_dimensions>
VoxelisedGeometricalInfo<num_dimensions>::
VoxelisedGeometricalInfo(
	const Offset& offset, const Spacing& spacing,
	const Size& size, const DirectionMatrix& direction)
	:
	offset(offset),
	spacing(spacing),
	size(size),
	direction(direction)
{}

template <int num_dimensions>
const typename VoxelisedGeometricalInfo<num_dimensions>::Offset
VoxelisedGeometricalInfo<num_dimensions>::get_offset() const
{
	return offset;
}

template <int num_dimensions>
const typename VoxelisedGeometricalInfo<num_dimensions>::Spacing
VoxelisedGeometricalInfo<num_dimensions>::get_spacing() const
{
	return spacing;
}

template <int num_dimensions>
const typename VoxelisedGeometricalInfo<num_dimensions>::Size
VoxelisedGeometricalInfo<num_dimensions>::get_size() const
{
	return size;
}

template <int num_dimensions>
const typename VoxelisedGeometricalInfo<num_dimensions>::DirectionMatrix
VoxelisedGeometricalInfo<num_dimensions>::get_direction() const
{
	return direction;
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
				direction[dim][axis] * spacing[dim];
		}
		index_to_physical_point_matrix[dim][num_dimensions] = offset[dim];
		index_to_physical_point_matrix[num_dimensions][dim] = 0;
	}
	index_to_physical_point_matrix[num_dimensions][num_dimensions] = 1;
	return index_to_physical_point_matrix;
}

}

#endif
