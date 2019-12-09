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

#include <array>

namespace sirf {

template <int num_physical_dimensions, int num_index_dimensions>
class GeometricalInfo {
public:
	typedef std::array<float, num_physical_dimensions>     Coordinate;
	typedef std::array<unsigned int, num_index_dimensions> Index;
	// Eventually something here like
	// Coordinate transform_index_to_physical_point(Index)
	// Index transform_physical_point_to_index(Coordinate)

    /// Print info
    virtual void print_info() const = 0;
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
	typedef std::array<Coordinate, num_dimensions> DirectionMatrix;
	/*!
	   Each vector in Direction tells the direction of the axis in LPS
	   physical space.
	*/
	typedef std::array<std::array<float, num_dimensions+1>, num_dimensions+1>
		TransformMatrix;


	VoxelisedGeometricalInfo(
		const Offset& _offset, const Spacing& _spacing,
		const Size& _size, const DirectionMatrix& _direction);
	virtual ~VoxelisedGeometricalInfo() {};

    const Offset get_offset() const;
    const Spacing get_spacing() const;
    const Size get_size() const;
    const DirectionMatrix get_direction() const;

    const TransformMatrix calculate_index_to_physical_point_matrix() const;

    /// Print info
    virtual void print_info() const;

private:
	Offset _offset;
	Spacing _spacing;
	Size _size;
	DirectionMatrix _direction;
};

typedef GeometricalInfo<3, 3> GeometricalInfo3D;
typedef VoxelisedGeometricalInfo<3> VoxelisedGeometricalInfo3D;
typedef VoxelisedGeometricalInfo<3>::TransformMatrix TransformMatrix3D;

}

#endif
