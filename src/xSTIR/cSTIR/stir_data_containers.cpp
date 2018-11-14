/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC

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

#include "stir_data_containers.h"

using namespace stir;
using namespace sirf;

std::string PETAcquisitionData::_storage_scheme;
shared_ptr<PETAcquisitionData> PETAcquisitionData::_template;

float
PETAcquisitionData::norm()
{
	double t = 0.0;
	for (int s = 0; s <= get_max_segment_num(); ++s)
	{
		SegmentBySinogram<float> seg = get_segment_by_sinogram(s);
		SegmentBySinogram<float>::full_iterator seg_iter;
		for (seg_iter = seg.begin_all(); seg_iter != seg.end_all();) {
			double r = *seg_iter++;
			t += r*r;
		}
		if (s != 0) {
			seg = get_segment_by_sinogram(-s);
			for (seg_iter = seg.begin_all(); seg_iter != seg.end_all();) {
				double r = *seg_iter++;
				t += r*r;
			}
		}
	}
	return sqrt((float)t);
}

//void
//PETAcquisitionData::mult(float a, const aDataContainer<float>& a_x)
//{
//	PETAcquisitionData& x = (PETAcquisitionData&)a_x;
//	int n = get_max_segment_num();
//	int nx = x.get_max_segment_num();
//	for (int s = 0; s <= n && s <= nx; ++s)
//	{
//		SegmentBySinogram<float> seg = get_empty_segment_by_sinogram(s);
//		SegmentBySinogram<float> sx = x.get_segment_by_sinogram(s);
//		SegmentBySinogram<float>::full_iterator seg_iter;
//		SegmentBySinogram<float>::full_iterator sx_iter;
//		for (seg_iter = seg.begin_all(), sx_iter = sx.begin_all();
//			seg_iter != seg.end_all() && sx_iter != sx.end_all();
//			/*empty*/) {
//			*seg_iter++ = float(a*double(*sx_iter++));
//		}
//		set_segment(seg);
//		if (s != 0) {
//			seg = get_empty_segment_by_sinogram(-s);
//			sx = x.get_segment_by_sinogram(-s);
//			for (seg_iter = seg.begin_all(), sx_iter = sx.begin_all();
//				seg_iter != seg.end_all() && sx_iter != sx.end_all();
//				/*empty*/)
//				*seg_iter++ = float(a*double(*sx_iter++));
//			set_segment(seg);
//		}
//	}
//}

float
PETAcquisitionData::dot(const aDataContainer<float>& a_x)
{
	PETAcquisitionData& x = (PETAcquisitionData&)a_x;
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
	return (float)t;
}

void
PETAcquisitionData::inv(float amin, const aDataContainer<float>& a_x)
{
	PETAcquisitionData& x = (PETAcquisitionData&)a_x;
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
PETAcquisitionData::axpby(
float a, const aDataContainer<float>& a_x,
float b, const aDataContainer<float>& a_y
)
{
	PETAcquisitionData& x = (PETAcquisitionData&)a_x;
	PETAcquisitionData& y = (PETAcquisitionData&)a_y;
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
PETAcquisitionData::multiply(
const aDataContainer<float>& a_x,
const aDataContainer<float>& a_y
)
{
	PETAcquisitionData& x = (PETAcquisitionData&)a_x;
	PETAcquisitionData& y = (PETAcquisitionData&)a_y;
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
const aDataContainer<float>& a_x,
const aDataContainer<float>& a_y
)
{
	PETAcquisitionData& x = (PETAcquisitionData&)a_x;
	PETAcquisitionData& y = (PETAcquisitionData&)a_y;
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

float
PETImageData::norm()
{
#ifdef _MSC_VER
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

float
PETImageData::dot(const aDataContainer<float>& a_x)
{
	PETImageData& x = (PETImageData&)a_x;
#ifdef _MSC_VER
	Image3DF::full_iterator iter;
	Image3DF::const_full_iterator iter_x;
#else
	typename Array<3, float>::full_iterator iter;
	typename Array<3, float>::const_full_iterator iter_x;
#endif

	double s = 0.0;
	for (iter = data().begin_all(), iter_x = x.data().begin_all();
		iter != data().end_all() && iter_x != x.data().end_all(); iter++, iter_x++) {
		double t = *iter;
		s += t * (*iter_x);
	}
	return (float)s;
}

//void
//PETImageData::mult(float a, const aDataContainer<float>& a_x)
//{
//	PETImageData& x = (PETImageData&)a_x;
//#ifdef _MSC_VER
//	Image3DF::full_iterator iter;
//	Image3DF::const_full_iterator iter_x;
//#else
//	typename Array<3, float>::full_iterator iter;
//	typename Array<3, float>::const_full_iterator iter_x;
//#endif
//
//	for (iter = data().begin_all(), iter_x = x.data().begin_all();
//		iter != data().end_all() && iter_x != x.data().end_all(); iter++, iter_x++)
//		*iter = a * (*iter_x);
//}

void
PETImageData::multiply(
const aDataContainer<float>& a_x,
const aDataContainer<float>& a_y)
{
	PETImageData& x = (PETImageData&)a_x;
	PETImageData& y = (PETImageData&)a_y;
#ifdef _MSC_VER
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
PETImageData::divide(
const aDataContainer<float>& a_x,
const aDataContainer<float>& a_y)
{
	PETImageData& x = (PETImageData&)a_x;
	PETImageData& y = (PETImageData&)a_y;
#ifdef _MSC_VER
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
		THROW("division by zero in PETImageData::divide");

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

void
PETImageData::axpby(
float a, const aDataContainer<float>& a_x,
float b, const aDataContainer<float>& a_y)
{
	PETImageData& x = (PETImageData&)a_x;
	PETImageData& y = (PETImageData&)a_y;
#ifdef _MSC_VER
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

int
PETImageData::get_dimensions(int* dim) const
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
PETImageData::get_voxel_sizes(float* vsize) const
{
	const Voxels3DF& voxels = (const Voxels3DF&)*_data;
	CartesianCoordinate3D<float> vs = voxels.get_voxel_size();
	for (int i = 0; i < 3; i++)
		vsize[i] = vs[i + 1];
}

int
PETImageData::get_data(float* data) const
{
	Image3DF& image = *_data;
	Coordinate3D<int> min_indices;
	Coordinate3D<int> max_indices;
	if (!image.get_regular_range(min_indices, max_indices))
		return -1;
	std::copy(image.begin_all(), image.end_all(), data);
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
	return 0;
}

int
PETImageData::set_data(const float* data)
{
	Image3DF& image = *_data;
	Coordinate3D<int> min_indices;
	Coordinate3D<int> max_indices;
	if (!image.get_regular_range(min_indices, max_indices))
		return -1;
	for (int z = min_indices[1], i = 0; z <= max_indices[1]; z++) {
		for (int y = min_indices[2]; y <= max_indices[2]; y++) {
			for (int x = min_indices[3]; x <= max_indices[3]; x++, i++) {
				image[z][y][x] = data[i];
			}
		}
	}
	return 0;
}
