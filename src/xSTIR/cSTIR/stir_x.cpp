/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
Copyright 2015 - 2017 University College London.

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

#include "stir_types.h"
#include "stir_x.h"

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

void
PETAcquisitionData::mult(float a, const aDataContainer<float>& a_x)
{
	PETAcquisitionData& x = (PETAcquisitionData&)a_x;
	int n = get_max_segment_num();
	int nx = x.get_max_segment_num();
	for (int s = 0; s <= n && s <= nx; ++s)
	{
		SegmentBySinogram<float> seg = get_empty_segment_by_sinogram(s);
		SegmentBySinogram<float> sx = x.get_segment_by_sinogram(s);
		SegmentBySinogram<float>::full_iterator seg_iter;
		SegmentBySinogram<float>::full_iterator sx_iter;
		for (seg_iter = seg.begin_all(), sx_iter = sx.begin_all();
			seg_iter != seg.end_all() && sx_iter != sx.end_all();
			/*empty*/) {
			*seg_iter++ = float(a*double(*sx_iter++));
		}
		set_segment(seg);
		if (s != 0) {
			seg = get_segment_by_sinogram(-s);
			sx = x.get_segment_by_sinogram(-s);
			for (seg_iter = seg.begin_all(), sx_iter = sx.begin_all();
				seg_iter != seg.end_all() && sx_iter != sx.end_all();
				/*empty*/)
				*seg_iter++ = float(a*double(*sx_iter++));
			set_segment(seg);
		}
	}
}

void
PETAcquisitionData::inv(float amin, const aDataContainer<float>& a_x)
{
	PETAcquisitionData& x = (PETAcquisitionData&)a_x;
	int n = get_max_segment_num();
	int nx = x.get_max_segment_num();
	for (int s = 0; s <= n && s <= nx; ++s)
	{
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
			seg = get_segment_by_sinogram(-s);
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
			seg = get_segment_by_sinogram(-s);
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
		iter != data().end_all(), iter_x != x.data().end_all(); iter++, iter_x++) {
		double t = *iter;
		s += t * (*iter_x);
	}
	return (float)s;
}

void 
PETImageData::mult(float a, const aDataContainer<float>& a_x)
{
	PETImageData& x = (PETImageData&)a_x;
#ifdef _MSC_VER
	Image3DF::full_iterator iter;
	Image3DF::const_full_iterator iter_x;
#else
	typename Array<3, float>::full_iterator iter;
	typename Array<3, float>::const_full_iterator iter_x;
#endif

	for (iter = data().begin_all(), iter_x = x.data().begin_all();
		iter != data().end_all(), iter_x != x.data().end_all(); iter++, iter_x++)
		*iter = a * (*iter_x);
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
		iter != data().end_all(), 
		iter_x != x.data().end_all(), iter_y != y.data().end_all(); 
		iter++, iter_x++, iter_y++)
		*iter = a * (*iter_x) + b * (*iter_y);
}
