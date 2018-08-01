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

template <int num_dimensions, typename CoordT>
class aCoordinate {
	typedef CoordT _CoordsT[num_dimensions];

public:
	inline CoordT& operator[](const int d)
		{ return coords[d]; }
	inline CoordT operator[](const int d) const
		{ return coords[d]; }

private:
	_CoordsT coords;
};


template <int num_physical_dimensions, int num_index_dimensions>
class GeometricalInfo {
	typedef aCoordinate<num_physical_dimensions, float>     Coordinate;
	typedef aCoordinate<num_index_dimensions, unsigned int> Index;
	// Eventually something here like
	// Coordinate transform_index_to_physical_point(Index)
	// Index transform_physical_point_to_index(Coordinate)
};


template <int num_dimensions>
class VoxelisedGeometricalInfo :
	public GeometricalInfo<num_dimensions, num_dimensions> {
public:
	typedef aCoordinate<num_dimensions, float>        Offset;
	typedef aCoordinate<num_dimensions, float>        Spacing;
	typedef aCoordinate<num_dimensions, unsigned int> Size;
	typedef aCoordinate<num_dimensions, aCoordinate<num_dimensions, float> >
		Direction;


	VoxelisedGeometricalInfo(
		const Offset& offset, const Spacing& spacing,
		const Size& size, const Direction& direction);
	// GeometricalInfo(const GeometricalInfo& other) {
	// 	*this = other;
	// }
	virtual ~VoxelisedGeometricalInfo() {};

	const Offset get_offset();
	const Spacing get_spacing();
	const Size get_size();
	const Direction get_direction();

private:
	Offset offset;
	Spacing spacing;
	Size size;
	Direction direction;
};

typedef VoxelisedGeometricalInfo<3> VoxelisedGeometricalInfo3D;

#endif
