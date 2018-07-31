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
class Coordinate {
	typedef CoordT _CoordsT[num_dimensions];

public:
	inline CoordT& operator[](const int d)
		{ return coords[d]; }
	inline CoordT operator[](const int d) const
		{ return coords[d]; }

private:
	_CoordsT coords;
};

class GeometricalInfo {
private:
	const static int DIM = 3;

public:
	typedef Coordinate<DIM, float>                   Offset;
	typedef Coordinate<DIM, float>                   Spacing;
	typedef Coordinate<DIM, int>                     Size;
	typedef Coordinate<DIM, Coordinate<DIM, float> > Direction;


	GeometricalInfo(
		const Offset& offset, const Spacing& spacing,
		const Size& size, const Direction& direction);
	// GeometricalInfo(const GeometricalInfo& other) {
	// 	*this = other;
	// }
	virtual ~GeometricalInfo() {};
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

#endif
