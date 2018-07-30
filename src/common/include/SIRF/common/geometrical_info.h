#ifndef SIRF_GEOMETRICAL_INFO_TYPE
#define SIRF_GEOMETRICAL_INFO_TYPE

const int DIM = 3;

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
public:
	typedef Coordinate<DIM, float>                  Offset;
	typedef Coordinate<DIM, float>                  Spacing;
	typedef Coordinate<DIM, int>                    Size;
	typedef Coordinate<DIM, Coordinate<DIM, float>> Direction;


	GeometricalInfo(Offset offset, Spacing spacing, Size size, Direction direction);
	GeometricalInfo(const GeometricalInfo& other) {
		*this = other;
	}
	~GeometricalInfo() {};
	Offset get_offset();
	Spacing get_spacing();
	Size get_size();
	Direction get_direction();
	//Transfrom GetTransfrom();

private:
	Offset offset;
	Spacing spacing;
	Size size;
	Direction direction;
	//Transform transform;
};

#endif
