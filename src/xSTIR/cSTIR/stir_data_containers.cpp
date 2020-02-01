/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
Copyright 2018 University College London

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

#include "sirf/STIR/stir_data_containers.h"
#include "stir/KeyParser.h"
#include "stir/is_null_ptr.h"
#include "stir/zoom.h"

using namespace stir;
using namespace sirf;

//#define DYNAMIC_CAST(T, X, Y) T& X = (T&)Y
#define DYNAMIC_CAST(T, X, Y) T& X = dynamic_cast<T&>(Y)

std::string PETAcquisitionData::_storage_scheme;
shared_ptr<PETAcquisitionData> PETAcquisitionData::_template;

float PETAcquisitionData::norm() const
{
	double t = 0.0;

	for (int s = 0; s <= get_max_segment_num(); ++s)
	{
		SegmentBySinogram<float> seg = get_segment_by_sinogram(s);
		SegmentBySinogram<float>::full_iterator seg_iter;

                for (seg_iter = seg.begin_all(); seg_iter != seg.end_all();)
                {
			double r = *seg_iter++;
			t += r*r;
		}

                if (s != 0)
                {
			seg = get_segment_by_sinogram(-s);

                        for (seg_iter = seg.begin_all(); seg_iter != seg.end_all();)
                        {
				double r = *seg_iter++;
				t += r*r;
			}
		}
	}

	return sqrt((float)t);
}

void
PETAcquisitionData::dot(const DataContainer& a_x, void* ptr) const
{
	//PETAcquisitionData& x = (PETAcquisitionData&)a_x;
	DYNAMIC_CAST(const PETAcquisitionData, x, a_x);
	int n = get_max_segment_num();
	int nx = x.get_max_segment_num();
	double t = 0;
	for (int s = 0; s <= n && s <= nx; ++s)
	{
		SegmentBySinogram<float> seg = get_segment_by_sinogram(s);
		SegmentBySinogram<float> sx = x.get_segment_by_sinogram(s);
		SegmentBySinogram<float>::full_iterator seg_iter;
		SegmentBySinogram<float>::full_iterator sx_iter;
		for (seg_iter = seg.begin_all(), sx_iter = sx.begin_all();
			seg_iter != seg.end_all() && sx_iter != sx.end_all();
			/*empty*/) {
			t += (*seg_iter++)*double(*sx_iter++);
		}
		if (s != 0) {
			seg = get_segment_by_sinogram(-s);
			sx = x.get_segment_by_sinogram(-s);
			for (seg_iter = seg.begin_all(), sx_iter = sx.begin_all();
				seg_iter != seg.end_all() && sx_iter != sx.end_all();
				/*empty*/)
				t += (*seg_iter++)*double(*sx_iter++);
		}
	}
	float* ptr_t = (float*)ptr;
	*ptr_t = (float)t;
}

void
PETAcquisitionData::axpby(
const void* ptr_a, const DataContainer& a_x,
const void* ptr_b, const DataContainer& a_y
)
{
	float a = *(float*)ptr_a;
	float b = *(float*)ptr_b;
	DYNAMIC_CAST(const PETAcquisitionData, x, a_x);
	DYNAMIC_CAST(const PETAcquisitionData, y, a_y);
	//PETAcquisitionData& x = (PETAcquisitionData&)a_x;
	//PETAcquisitionData& y = (PETAcquisitionData&)a_y;
	int n = get_max_segment_num();
	int nx = x.get_max_segment_num();
	int ny = y.get_max_segment_num();
	for (int s = 0; s <= n && s <= nx && s <= ny; ++s)
	{
		SegmentBySinogram<float> seg = get_empty_segment_by_sinogram(s);
		SegmentBySinogram<float> sx = x.get_segment_by_sinogram(s);
		SegmentBySinogram<float> sy = y.get_segment_by_sinogram(s);
		SegmentBySinogram<float>::full_iterator seg_iter;
		SegmentBySinogram<float>::full_iterator sx_iter;
		SegmentBySinogram<float>::full_iterator sy_iter;
		for (seg_iter = seg.begin_all(),
			sx_iter = sx.begin_all(), sy_iter = sy.begin_all();
			seg_iter != seg.end_all() &&
			sx_iter != sx.end_all() && sy_iter != sy.end_all();
		/*empty*/) {
			*seg_iter++ = float(a*double(*sx_iter++) + b*double(*sy_iter++));
		}
		set_segment(seg);
		if (s != 0) {
			seg = get_empty_segment_by_sinogram(-s);
			sx = x.get_segment_by_sinogram(-s);
			sy = y.get_segment_by_sinogram(-s);
			for (seg_iter = seg.begin_all(),
				sx_iter = sx.begin_all(), sy_iter = sy.begin_all();
				seg_iter != seg.end_all() &&
				sx_iter != sx.end_all() && sy_iter != sy.end_all();
			/*empty*/) {
				*seg_iter++ = float(a*double(*sx_iter++) + b*double(*sy_iter++));
			}
			set_segment(seg);
		}
	}
}

void
PETAcquisitionData::inv(float amin, const DataContainer& a_x)
{
	//PETAcquisitionData& x = (PETAcquisitionData&)a_x;
	DYNAMIC_CAST(const PETAcquisitionData, x, a_x);
	int n = get_max_segment_num();
	int nx = x.get_max_segment_num();
	for (int s = 0; s <= n && s <= nx; ++s)
	{
		//std::cout << "processing segment " << s << std::endl;
		SegmentBySinogram<float> seg = get_empty_segment_by_sinogram(s);
		SegmentBySinogram<float> sx = x.get_segment_by_sinogram(s);
		SegmentBySinogram<float>::full_iterator seg_iter;
		SegmentBySinogram<float>::full_iterator sx_iter;
		for (seg_iter = seg.begin_all(), sx_iter = sx.begin_all();
			seg_iter != seg.end_all() && sx_iter != sx.end_all();
			/*empty*/)
			*seg_iter++ = float(1.0 / std::max(amin, *sx_iter++));
		set_segment(seg);
		if (s != 0) {
			//std::cout << "processing segment " << -s << std::endl;
			seg = get_empty_segment_by_sinogram(-s);
			sx = x.get_segment_by_sinogram(-s);
			for (seg_iter = seg.begin_all(), sx_iter = sx.begin_all();
				seg_iter != seg.end_all() && sx_iter != sx.end_all();
				/*empty*/) {
				*seg_iter++ = float(1.0 / std::max(amin, *sx_iter++));
			}
			set_segment(seg);
		}
	}
}

void
PETAcquisitionData::multiply(
const DataContainer& a_x,
const DataContainer& a_y
)
{
	//PETAcquisitionData& x = (PETAcquisitionData&)a_x;
	//PETAcquisitionData& y = (PETAcquisitionData&)a_y;
	DYNAMIC_CAST(const PETAcquisitionData, x, a_x);
	DYNAMIC_CAST(const PETAcquisitionData, y, a_y);
	int n = get_max_segment_num();
	int nx = x.get_max_segment_num();
	int ny = y.get_max_segment_num();
	for (int s = 0; s <= n && s <= nx && s <= ny; ++s)
	{
		SegmentBySinogram<float> seg = get_empty_segment_by_sinogram(s);
		SegmentBySinogram<float> sx = x.get_segment_by_sinogram(s);
		SegmentBySinogram<float> sy = y.get_segment_by_sinogram(s);
		SegmentBySinogram<float>::full_iterator seg_iter;
		SegmentBySinogram<float>::full_iterator sx_iter;
		SegmentBySinogram<float>::full_iterator sy_iter;
		for (seg_iter = seg.begin_all(),
			sx_iter = sx.begin_all(), sy_iter = sy.begin_all();
			seg_iter != seg.end_all() &&
			sx_iter != sx.end_all() && sy_iter != sy.end_all();
		/*empty*/) {
			*seg_iter++ = (*sx_iter++) * (*sy_iter++);
		}
		set_segment(seg);
		if (s != 0) {
			seg = get_empty_segment_by_sinogram(-s);
			sx = x.get_segment_by_sinogram(-s);
			sy = y.get_segment_by_sinogram(-s);
			for (seg_iter = seg.begin_all(),
				sx_iter = sx.begin_all(), sy_iter = sy.begin_all();
				seg_iter != seg.end_all() &&
				sx_iter != sx.end_all() && sy_iter != sy.end_all();
			/*empty*/) {
				*seg_iter++ = (*sx_iter++) * (*sy_iter++);
			}
			set_segment(seg);
		}
	}
}

void
PETAcquisitionData::divide(
const DataContainer& a_x,
const DataContainer& a_y
)
{
	//PETAcquisitionData& x = (PETAcquisitionData&)a_x;
	//PETAcquisitionData& y = (PETAcquisitionData&)a_y;
	DYNAMIC_CAST(const PETAcquisitionData, x, a_x);
	DYNAMIC_CAST(const PETAcquisitionData, y, a_y);
	int n = get_max_segment_num();
	int nx = x.get_max_segment_num();
	int ny = y.get_max_segment_num();
	for (int s = 0; s <= n && s <= nx && s <= ny; ++s)
	{
		SegmentBySinogram<float> seg = get_empty_segment_by_sinogram(s);
		SegmentBySinogram<float> sx = x.get_segment_by_sinogram(s);
		SegmentBySinogram<float> sy = y.get_segment_by_sinogram(s);
		SegmentBySinogram<float>::full_iterator seg_iter;
		SegmentBySinogram<float>::full_iterator sx_iter;
		SegmentBySinogram<float>::full_iterator sy_iter;
		for (seg_iter = seg.begin_all(),
			sx_iter = sx.begin_all(), sy_iter = sy.begin_all();
			seg_iter != seg.end_all() &&
			sx_iter != sx.end_all() && sy_iter != sy.end_all();
		/*empty*/) {
			*seg_iter++ = (*sx_iter++) / (*sy_iter++);
		}
		set_segment(seg);
		if (s != 0) {
			seg = get_empty_segment_by_sinogram(-s);
			sx = x.get_segment_by_sinogram(-s);
			sy = y.get_segment_by_sinogram(-s);
			for (seg_iter = seg.begin_all(),
				sx_iter = sx.begin_all(), sy_iter = sy.begin_all();
				seg_iter != seg.end_all() &&
				sx_iter != sx.end_all() && sy_iter != sy.end_all();
			/*empty*/) {
				*seg_iter++ = (*sx_iter++) / (*sy_iter++);
			}
			set_segment(seg);
		}
	}
}

STIRImageData::STIRImageData(const ImageData& id)
{
    throw std::runtime_error("TODO - create STIRImageData from general SIRFImageData.");
    /* The following is incorrect.
    Dimensions dim = id.dimensions();
    int nx = dim["x"];
    int ny = dim["y"];
    int nz = 1;
    Dimensions::iterator it = dim.begin();
    while (it != dim.end()) {
        if (it->first != "x" && it->first != "y")
            nz *= it->second;
        ++it;
    }
    Voxels3DF voxels(stir::IndexRange3D(0, nz - 1,
        -(ny / 2), -(ny / 2) + ny - 1, -(nx / 2), -(nx / 2) + nx - 1),
        Coord3DF(0, 0, 0),
        Coord3DF(3, 3, 3.375));
    _data.reset(voxels.clone());
    copy(id.begin(), begin(), end());

    // Set up the geom info
    this->set_up_geom_info();*/
}

void
STIRImageData::write(const std::string &filename) const
{
    this->write(filename,"");
}

void
STIRImageData::write(const std::string &filename, const std::string &format_file) const
{
    const Image3DF& image = this->data();
    shared_ptr<OutputFileFormat<Image3DF> > format_sptr;

    if (!format_file.empty()) {
        KeyParser parser;
        parser.add_start_key("OutputFileFormat Parameters");
        parser.add_parsing_key("output file format type", &format_sptr);
        parser.add_stop_key("END");
        parser.parse(format_file.c_str());
        if(is_null_ptr(format_sptr))
            throw std::runtime_error("STIRImageData::write: Parsing of output format file (" + format_file + ") failed "
                                     "(see examples/parameter_files/STIR_output_file_format_xxx.par for help).");
    }
    if(is_null_ptr(format_sptr))
        format_sptr = OutputFileFormat<Image3DF>::default_sptr();

    format_sptr->write_to_file(filename, image);
}

void
STIRImageData::dot(const DataContainer& a_x, void* ptr) const
{
	//STIRImageData& x = (STIRImageData&)a_x;
	DYNAMIC_CAST(const STIRImageData, x, a_x);
#if defined(_MSC_VER) && _MSC_VER < 1900
	Image3DF::const_full_iterator iter;
	Image3DF::const_full_iterator iter_x;
#else
	typename Array<3, float>::const_full_iterator iter;
	typename Array<3, float>::const_full_iterator iter_x;
#endif

	double s = 0.0;
	for (iter = data().begin_all(), iter_x = x.data().begin_all();
		iter != data().end_all() && iter_x != x.data().end_all();
		iter++, iter_x++) {
		double t = *iter;
		s += t * (*iter_x);
	}
	float* ptr_s = (float*)ptr;
	*ptr_s = (float)s;
}

void
STIRImageData::axpby(
const void* ptr_a, const DataContainer& a_x,
const void* ptr_b, const DataContainer& a_y)
{
	float a = *(float*)ptr_a;
	float b = *(float*)ptr_b;
	DYNAMIC_CAST(const STIRImageData, x, a_x);
	DYNAMIC_CAST(const STIRImageData, y, a_y);
	//STIRImageData& x = (STIRImageData&)a_x;
	//STIRImageData& y = (STIRImageData&)a_y;
#if defined(_MSC_VER) && _MSC_VER < 1900
	Image3DF::full_iterator iter;
	Image3DF::const_full_iterator iter_x;
	Image3DF::const_full_iterator iter_y;
#else
	typename Array<3, float>::full_iterator iter;
	typename Array<3, float>::const_full_iterator iter_x;
	typename Array<3, float>::const_full_iterator iter_y;
#endif

	for (iter = data().begin_all(),
		iter_x = x.data().begin_all(), iter_y = y.data().begin_all();
		iter != data().end_all() &&
		iter_x != x.data().end_all() && iter_y != y.data().end_all();
	iter++, iter_x++, iter_y++)
		*iter = a * (*iter_x) + b * (*iter_y);
}

float
STIRImageData::norm() const
{
#if defined(_MSC_VER) && _MSC_VER < 1900
	//Array<3, float>::const_full_iterator iter;
	Image3DF::const_full_iterator iter;
#else
	typename Array<3, float>::const_full_iterator iter;
#endif
	double s = 0.0;
	int i = 0;
	for (iter = _data->begin_all(); iter != _data->end_all(); iter++, i++) {
		double t = *iter;
		s += t*t;
	}
	//std::cout << "voxels count: " << i << std::endl;
	return (float)sqrt(s);
}

void
STIRImageData::multiply(
const DataContainer& a_x,
const DataContainer& a_y)
{
	//STIRImageData& x = (STIRImageData&)a_x;
	//STIRImageData& y = (STIRImageData&)a_y;
	DYNAMIC_CAST(const STIRImageData, x, a_x);
	DYNAMIC_CAST(const STIRImageData, y, a_y);
#if defined(_MSC_VER) && _MSC_VER < 1900
	Image3DF::full_iterator iter;
	Image3DF::const_full_iterator iter_x;
	Image3DF::const_full_iterator iter_y;
#else
	typename Array<3, float>::full_iterator iter;
	typename Array<3, float>::const_full_iterator iter_x;
	typename Array<3, float>::const_full_iterator iter_y;
#endif

	for (iter = data().begin_all(),
		iter_x = x.data().begin_all(), iter_y = y.data().begin_all();
		iter != data().end_all() &&
		iter_x != x.data().end_all() && iter_y != y.data().end_all();
	iter++, iter_x++, iter_y++)
		*iter = (*iter_x) * (*iter_y);
}

void
STIRImageData::divide(
const DataContainer& a_x,
const DataContainer& a_y)
{
	//STIRImageData& x = (STIRImageData&)a_x;
	//STIRImageData& y = (STIRImageData&)a_y;
	DYNAMIC_CAST(const STIRImageData, x, a_x);
	DYNAMIC_CAST(const STIRImageData, y, a_y);
#if defined(_MSC_VER) && _MSC_VER < 1900
	Image3DF::full_iterator iter;
	Image3DF::const_full_iterator iter_x;
	Image3DF::const_full_iterator iter_y;
#else
	typename Array<3, float>::full_iterator iter;
	typename Array<3, float>::const_full_iterator iter_x;
	typename Array<3, float>::const_full_iterator iter_y;
#endif

	float vmax = 0.0;
	for (
		iter_x = x.data().begin_all(), iter_y = y.data().begin_all();
		iter_x != x.data().end_all() && iter_y != y.data().end_all();
	iter_x++, iter_y++) {
		float vy = abs(*iter_y);
		if (vy > vmax)
			vmax = vy;
	}
	float vmin = 1e-6*vmax;
	if (vmin == 0.0)
		THROW("division by zero in STIRImageData::divide");

	for (iter = data().begin_all(),
		iter_x = x.data().begin_all(), iter_y = y.data().begin_all();
		iter != data().end_all() &&
		iter_x != x.data().end_all() && iter_y != y.data().end_all();
	iter++, iter_x++, iter_y++) {
		float vy = *iter_y;
		if (vy >= 0 && vy < vmin)
			vy = vmin;
		else if (vy < 0 && vy > -vmin)
			vy = -vmin;
		*iter = (*iter_x) / vy;
	}
}

int
STIRImageData::get_dimensions(int* dim) const
{
	const Image3DF& image = *_data;
	dim[0] = 0;
	dim[1] = 0;
	dim[2] = 0;
	Coordinate3D<int> min_indices;
	Coordinate3D<int> max_indices;
	if (!image.get_regular_range(min_indices, max_indices))
		return -1;
	for (int i = 0; i < 3; i++)
		dim[i] = max_indices[i + 1] - min_indices[i + 1] + 1;
	return 0;
}

void
STIRImageData::get_voxel_sizes(float* vsize) const
{
	//const Voxels3DF& voxels = (const Voxels3DF&)*_data;
	DYNAMIC_CAST(const Voxels3DF, voxels, *_data);
	CartesianCoordinate3D<float> vs = voxels.get_voxel_size();
	for (int i = 0; i < 3; i++)
		vsize[i] = vs[i + 1];
}

void
STIRImageData::get_data(float* data) const
{
	Image3DF& image = *_data;
	Coordinate3D<int> min_indices;
	Coordinate3D<int> max_indices;
	if (!image.get_regular_range(min_indices, max_indices))
		throw LocalisedException("irregular STIR image", __FILE__, __LINE__);
		//return -1;
	//std::cout << "trying new const iterator...\n";
	STIRImageData::Iterator_const iter(begin());
	for (int i = 0; iter != end(); ++i, ++iter)
		data[i] = *iter;
	//std::copy(begin(), end(), data);
	//std::copy(image.begin_all(), image.end_all(), data);
	//auto iter = image.begin_all();
	//for (int i = 0; iter != image.end_all(); i++, iter++)
	//	data[i] = *iter;
	//for (int z = min_indices[1], i = 0; z <= max_indices[1]; z++) {
	//	for (int y = min_indices[2]; y <= max_indices[2]; y++) {
	//		for (int x = min_indices[3]; x <= max_indices[3]; x++, i++) {
	//			data[i] = image[z][y][x];
	//		}
	//	}
	//}
	//return 0;
}

void
STIRImageData::set_data(const float* data)
{
	Image3DF& image = *_data;
	Coordinate3D<int> min_indices;
	Coordinate3D<int> max_indices;
	if (!image.get_regular_range(min_indices, max_indices))
		throw LocalisedException("irregular STIR image", __FILE__, __LINE__);
	//return -1;
	size_t n = 1;
	for (int i = 0; i < 3; i++)
		n *= (max_indices[i + 1] - min_indices[i + 1] + 1);
	//std::cout << "trying new iterator...\n";
	STIRImageData::Iterator iter(begin());
	for (int i = 0; iter != end(); ++i, ++iter)
		*iter = data[i];
	//std::copy(data, data + n, begin());
	//std::copy(data, data + n, image.begin_all());
	//for (int z = min_indices[1], i = 0; z <= max_indices[1]; z++) {
	//	for (int y = min_indices[2]; y <= max_indices[2]; y++) {
	//		for (int x = min_indices[3]; x <= max_indices[3]; x++, i++) {
	//			image[z][y][x] = data[i];
	//		}
	//	}
	//}
	//return 0;
}

void
STIRImageData::
zoom_image(const Coord3DF &zooms, const Coord3DF &offsets_in_mm,
           const Coord3DI &new_sizes, const char *zoom_options_str)
{
    stir::ZoomOptions zoom_options;
    if (strcmp(zoom_options_str,"preserve_sum")==0)
        zoom_options = stir::ZoomOptions::preserve_sum;
    else if (strcmp(zoom_options_str,"preserve_values")==0)
        zoom_options = stir::ZoomOptions::preserve_values;
    else if (strcmp(zoom_options_str,"preserve_projections")==0)
        zoom_options = stir::ZoomOptions::preserve_projections;
    else
        throw std::runtime_error("zoom_image: unknown scaling option - " + std::string(zoom_options_str));

    this->zoom_image(zooms, offsets_in_mm, new_sizes, zoom_options);
}

void
STIRImageData::
zoom_image(const Coord3DF &zooms, const Coord3DF &offsets_in_mm,
           const Coord3DI &new_sizes_in, const stir::ZoomOptions zoom_options)
{
    // We need the underyling image as a VoxelsOnCartesianGrid
    DYNAMIC_CAST(Voxels3DF, voxels, this->data());

    int dim[3];
    this->get_dimensions(dim);

    // If any sizes have been set to <= 0, set to image size
    Coord3DI new_sizes(new_sizes_in);
    for (unsigned i=0; i<3; ++i)
        if (new_sizes[int(i)]<=0)
            new_sizes[int(i)] = dim[i];

    // Zoom the image
    voxels = stir::zoom_image(voxels, zooms, offsets_in_mm, new_sizes, zoom_options);
}

void
STIRImageData::
move_to_scanner_centre(const PETAcquisitionData &)
{
    this->_data->set_origin(CartesianCoordinate3D<float>{0.f,0.f,0.f});
}

void
STIRImageData::set_up_geom_info()
{
    const Voxels3DF* const vox_image = dynamic_cast<const Voxels3DF*>(&data());

    // If cast failed, throw error
    if (!vox_image)
        throw std::runtime_error("Can't determine geometry for this image type");

    // SIRF offest is STIR's LPS location of first voxel
    VoxelisedGeometricalInfo3D::Offset offset;
    const stir::CartesianCoordinate3D<int> first_vox
        = vox_image->get_min_indices();
    const Coord3DF first_vox_coord
        = vox_image->get_LPS_coordinates_for_indices(first_vox);
    offset[0] = first_vox_coord.x();
    offset[1] = first_vox_coord.y();
    offset[2] = first_vox_coord.z();

    // SIRF and STIR share size definition
    VoxelisedGeometricalInfo3D::Size size;
    size[0] = vox_image->get_x_size();
    size[1] = vox_image->get_y_size();
    size[2] = vox_image->get_z_size();

    // SIRF's spacing is STIR's voxel size, but with different order
    VoxelisedGeometricalInfo3D::Spacing spacing;
    spacing[0] = vox_image->get_voxel_size()[3];
    spacing[1] = vox_image->get_voxel_size()[2];
    spacing[2] = vox_image->get_voxel_size()[1];

    // Find axes direction as the normalised vector between voxels
    // in each direction
    VoxelisedGeometricalInfo3D::DirectionMatrix direction;
    for (int axis = 0; axis < 3; axis++) {
        Coord3DI next_vox_along_axis(first_vox);
        next_vox_along_axis[3 - axis] += 1;
        const Coord3DF next_vox_along_axis_coord
            = vox_image->get_LPS_coordinates_for_indices(next_vox_along_axis);
        Coord3DF axis_direction
            = next_vox_along_axis_coord - first_vox_coord;
        axis_direction /= stir::norm(axis_direction);
        for (int dim = 0; dim < 3; dim++)
            direction[dim][axis] = axis_direction[3 - dim];
    }

    // Initialise the geom info shared pointer
    _geom_info_sptr = std::make_shared<VoxelisedGeometricalInfo3D>
                (offset,spacing,size,direction);
}
