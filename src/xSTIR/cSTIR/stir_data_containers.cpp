/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2017 - 2023 Rutherford Appleton Laboratory STFC
Copyright 2018 - 2023 University College London

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

#include "sirf/STIR/stir_data_containers.h"
#include "stir/KeyParser.h"
#include "stir/is_null_ptr.h"
#include "stir/zoom.h"
#include "stir/CartesianCoordinate3D.h"
#include "stir/numerics/norm.h"

using namespace stir;
using namespace sirf;

#define SIRF_DYNAMIC_CAST(T, X, Y) T& X = dynamic_cast<T&>(Y)

std::string STIRAcquisitionData::_storage_scheme;
std::shared_ptr<STIRAcquisitionData> STIRAcquisitionData::_template;

#ifdef STIR_TOF
#define TOF_LOOP  for (int k=data()->get_min_tof_pos_num(); k<=data()->get_max_tof_pos_num(); ++k)
#define TOF_ARG , k
#else
#define TOF_LOOP
#define TOF_ARG
#endif

float
STIRAcquisitionData::norm() const
{
#if STIR_VERSION >= 060200
        return data()->norm();
#else
	double t = 0.0;
        TOF_LOOP
          for (int s = data()->get_min_segment_num(); s <= data()->get_max_segment_num(); ++s)
	{
		const auto seg = get_segment_by_sinogram(s TOF_ARG);
		t += stir::norm_squared(seg.begin_all(), seg.end_all());
	}
	return static_cast<float>(std::sqrt(t));
#endif
}

void
STIRAcquisitionData::sum(void* ptr) const
{
	float* ptr_t = static_cast<float*>(ptr);
#if STIR_VERSION >= 060200
        *ptr_t = data()->sum();
#else
	int n = get_max_segment_num();
	double t = 0;
        TOF_LOOP
	for (int s = 0; s <= n; ++s)
	{
		SegmentBySinogram<float> seg = get_segment_by_sinogram(s TOF_ARG);
		SegmentBySinogram<float>::full_iterator seg_iter;
		for (seg_iter = seg.begin_all(); seg_iter != seg.end_all();)
			t += *seg_iter++;
		if (s != 0) {
			seg = get_segment_by_sinogram(-s TOF_ARG);
			for (seg_iter = seg.begin_all(); seg_iter != seg.end_all();)
				t += *seg_iter++;
		}
	}
	*ptr_t = (float)t;
#endif
}

void
STIRAcquisitionData::max(void* ptr) const
{
	float* ptr_t = static_cast<float*>(ptr);
#if STIR_VERSION >= 060200
        *ptr_t = data()->find_max();
#else
	int n = get_max_segment_num();
	float t = 0;
	bool init = true;
        TOF_LOOP
	for (int s = 0; s <= n; ++s)
	{
		SegmentBySinogram<float> seg = get_segment_by_sinogram(s TOF_ARG);
		SegmentBySinogram<float>::full_iterator seg_iter;
		for (seg_iter = seg.begin_all(); seg_iter != seg.end_all();)
			if (init) {
				t = *seg_iter++;
				init = false;
			}
			else
				t = std::max(t, *seg_iter++);
		if (s != 0) {
			seg = get_segment_by_sinogram(-s TOF_ARG);
			for (seg_iter = seg.begin_all(); seg_iter != seg.end_all();)
				t = std::max(t, *seg_iter++);
		}
	}
	*ptr_t = (float)t;
#endif
}

void
STIRAcquisitionData::min(void* ptr) const
{
	float* ptr_t = static_cast<float*>(ptr);
#if STIR_VERSION >= 060200
        *ptr_t = data()->find_min();
#else
	float t = 0;
	bool init = true;
        TOF_LOOP
        for (int s = data()->get_min_segment_num(); s <= get_max_segment_num(); ++s)
	{
                const auto t_seg = get_segment_by_sinogram(s TOF_ARG).find_min();
                if (init) {
                        init = false;
                        t = t_seg;
                }
                else {
                        t = std::min(t, t_seg);
		}
	}
	*ptr_t = (float)t;
#endif
}

void
STIRAcquisitionData::dot(const DataContainer& a_x, void* ptr) const
{
	SIRF_DYNAMIC_CAST(const STIRAcquisitionData, x, a_x);
	int n = get_max_segment_num();
	int nx = x.get_max_segment_num();
	double t = 0;
        TOF_LOOP
	for (int s = 0; s <= n && s <= nx; ++s)
	{
		SegmentBySinogram<float> seg = get_segment_by_sinogram(s TOF_ARG);
		SegmentBySinogram<float> sx = x.get_segment_by_sinogram(s TOF_ARG);
		SegmentBySinogram<float>::full_iterator seg_iter;
		SegmentBySinogram<float>::full_iterator sx_iter;
		for (seg_iter = seg.begin_all(), sx_iter = sx.begin_all();
			seg_iter != seg.end_all() && sx_iter != sx.end_all();
			/*empty*/) {
			t += (*seg_iter++) * double(*sx_iter++);
		}
		if (s != 0) {
			seg = get_segment_by_sinogram(-s TOF_ARG);
			sx = x.get_segment_by_sinogram(-s TOF_ARG);
			for (seg_iter = seg.begin_all(), sx_iter = sx.begin_all();
				seg_iter != seg.end_all() && sx_iter != sx.end_all();
				/*empty*/)
				t += (*seg_iter++) * double(*sx_iter++);
		}
	}
	float* ptr_t = static_cast<float*>(ptr);
	*ptr_t = (float)t;
}

void
STIRAcquisitionData::axpby(
const void* ptr_a, const DataContainer& a_x,
const void* ptr_b, const DataContainer& a_y
)
{
	//Add deprecation warning
    STIRAcquisitionData::xapyb(a_x, ptr_a, a_y, ptr_b);
}

void
STIRAcquisitionData::xapyb(
const DataContainer& a_x, const void* ptr_a,
const DataContainer& a_y, const void* ptr_b
)
{
    // Cast to correct types
	float a = *static_cast<const float*>(ptr_a);
	float b = *static_cast<const float*>(ptr_b);
	auto x = dynamic_cast<const STIRAcquisitionData*>(&a_x);
    auto y = dynamic_cast<const STIRAcquisitionData*>(&a_y);

    if (is_null_ptr(x) || is_null_ptr(x->data()) ||
            is_null_ptr(y) || is_null_ptr(y->data()))
        throw std::runtime_error("STIRAcquisitionData::xapyb: At least one argument is not"
                                 "STIRAcquisitionData or is not initialised.");

    // Call STIR's xapyb
    data()->xapyb(*x->data(), a, *y->data(), b);
}

void
STIRAcquisitionData::xapyb(
const DataContainer& a_x, const DataContainer& a_a,
const DataContainer& a_y, const DataContainer& a_b
)
{
    // Cast to correct types
    auto a = dynamic_cast<const STIRAcquisitionData*>(&a_a);
    auto b = dynamic_cast<const STIRAcquisitionData*>(&a_b);
    auto x = dynamic_cast<const STIRAcquisitionData*>(&a_x);
    auto y = dynamic_cast<const STIRAcquisitionData*>(&a_y);

    if (is_null_ptr(x) || is_null_ptr(x->data()) ||
            is_null_ptr(y) || is_null_ptr(y->data()) ||
            is_null_ptr(a) || is_null_ptr(a->data()) ||
            is_null_ptr(b) || is_null_ptr(b->data()))
        throw std::runtime_error("STIRAcquisitionData::xapyb: At least one argument is not"
                                 "STIRAcquisitionData or is not initialised.");

    // Call STIR's xapyb
    data()->xapyb(*x->data(), *a->data(), *y->data(), *b->data());
}

void
STIRAcquisitionData::inv(float amin, const DataContainer& a_x)
{
	SIRF_DYNAMIC_CAST(const STIRAcquisitionData, x, a_x);
	int n = get_max_segment_num();
	int nx = x.get_max_segment_num();
        TOF_LOOP
	for (int s = 0; s <= n && s <= nx; ++s)
	{
		//std::cout << "processing segment " << s << std::endl;
		SegmentBySinogram<float> seg = get_empty_segment_by_sinogram(s TOF_ARG);
		SegmentBySinogram<float> sx = x.get_segment_by_sinogram(s TOF_ARG);
		SegmentBySinogram<float>::full_iterator seg_iter;
		SegmentBySinogram<float>::full_iterator sx_iter;
		for (seg_iter = seg.begin_all(), sx_iter = sx.begin_all();
			seg_iter != seg.end_all() && sx_iter != sx.end_all();
			/*empty*/)
			*seg_iter++ = float(1.0 / std::max(amin, *sx_iter++));
		set_segment(seg);
		if (s != 0) {
			//std::cout << "processing segment " << -s << std::endl;
			seg = get_empty_segment_by_sinogram(-s TOF_ARG);
			sx = x.get_segment_by_sinogram(-s TOF_ARG);
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
STIRAcquisitionData::unary_op(
	const DataContainer& a_x,
	float(*f)(float)
)
{
	SIRF_DYNAMIC_CAST(const STIRAcquisitionData, x, a_x);
	int n = get_max_segment_num();
	int nx = x.get_max_segment_num();
        TOF_LOOP
	for (int s = 0; s <= n && s <= nx; ++s) {
		SegmentBySinogram<float> seg = get_empty_segment_by_sinogram(s TOF_ARG);
		SegmentBySinogram<float> sx = x.get_segment_by_sinogram(s TOF_ARG);
		SegmentBySinogram<float>::full_iterator seg_iter;
		SegmentBySinogram<float>::full_iterator sx_iter;
		for (seg_iter = seg.begin_all(), sx_iter = sx.begin_all();
			seg_iter != seg.end_all() && sx_iter != sx.end_all(); /*empty*/)
			*seg_iter++ = f(*sx_iter++);
		set_segment(seg);
		if (s > 0) {
			SegmentBySinogram<float> seg = get_empty_segment_by_sinogram(-s TOF_ARG);
			SegmentBySinogram<float> sx = x.get_segment_by_sinogram(-s TOF_ARG);
			SegmentBySinogram<float>::full_iterator seg_iter;
			SegmentBySinogram<float>::full_iterator sx_iter;
			for (seg_iter = seg.begin_all(), sx_iter = sx.begin_all();
				seg_iter != seg.end_all() && sx_iter != sx.end_all(); /*empty*/)
				*seg_iter++ = f(*sx_iter++);
			set_segment(seg);
		}
	}

}

void
STIRAcquisitionData::semibinary_op(
	const DataContainer& a_x,
	float y,
	float(*f)(float, float)
)
{
	SIRF_DYNAMIC_CAST(const STIRAcquisitionData, x, a_x);
	int n = get_max_segment_num();
	int nx = x.get_max_segment_num();
        TOF_LOOP
	for (int s = 0; s <= n && s <= nx; ++s) {
		SegmentBySinogram<float> seg = get_empty_segment_by_sinogram(s TOF_ARG);
		SegmentBySinogram<float> sx = x.get_segment_by_sinogram(s TOF_ARG);
		SegmentBySinogram<float>::full_iterator seg_iter;
		SegmentBySinogram<float>::full_iterator sx_iter;
		for (seg_iter = seg.begin_all(), sx_iter = sx.begin_all();
			seg_iter != seg.end_all() && sx_iter != sx.end_all(); /*empty*/)
			*seg_iter++ = f(*sx_iter++, y);
		set_segment(seg);
		if (s > 0) {
			SegmentBySinogram<float> seg = get_empty_segment_by_sinogram(-s TOF_ARG);
			SegmentBySinogram<float> sx = x.get_segment_by_sinogram(-s TOF_ARG);
			SegmentBySinogram<float>::full_iterator seg_iter;
			SegmentBySinogram<float>::full_iterator sx_iter;
			for (seg_iter = seg.begin_all(), sx_iter = sx.begin_all();
				seg_iter != seg.end_all() && sx_iter != sx.end_all(); /*empty*/)
				*seg_iter++ = f(*sx_iter++, y);
			set_segment(seg);
		}
	}

}

void
STIRAcquisitionData::binary_op(
	const DataContainer& a_x,
	const DataContainer& a_y,
	float(*f)(float, float)
)
{
	SIRF_DYNAMIC_CAST(const STIRAcquisitionData, x, a_x);
	SIRF_DYNAMIC_CAST(const STIRAcquisitionData, y, a_y);
	int n = get_max_segment_num();
	int nx = x.get_max_segment_num();
	int ny = y.get_max_segment_num();
	if (n != nx || n != ny)
		throw std::runtime_error("binary_op error: operands sizes differ");
	SegmentBySinogram<float>::full_iterator seg_iter;
	SegmentBySinogram<float>::full_iterator sx_iter;
	SegmentBySinogram<float>::full_iterator sy_iter;
        TOF_LOOP
	for (int s = 0; s <= n; ++s)
	{
		SegmentBySinogram<float> seg = get_empty_segment_by_sinogram(s TOF_ARG);
		SegmentBySinogram<float> sx = x.get_segment_by_sinogram(s TOF_ARG);
		SegmentBySinogram<float> sy = y.get_segment_by_sinogram(s TOF_ARG);
		if (seg.size_all() != sx.size_all() || seg.size_all() != sy.size_all())
			throw std::runtime_error("binary_op error: operands sizes differ");
		for (seg_iter = seg.begin_all(),
			sx_iter = sx.begin_all(), sy_iter = sy.begin_all();
			seg_iter != seg.end_all(); /*empty*/)
			*seg_iter++ = f(*sx_iter++, *sy_iter++);
		set_segment(seg);
		if (s != 0) {
			seg = get_empty_segment_by_sinogram(-s TOF_ARG);
			sx = x.get_segment_by_sinogram(-s TOF_ARG);
			sy = y.get_segment_by_sinogram(-s TOF_ARG);
			for (seg_iter = seg.begin_all(),
				sx_iter = sx.begin_all(), sy_iter = sy.begin_all();
				seg_iter != seg.end_all();	/*empty*/)
				*seg_iter++ = f(*sx_iter++, *sy_iter++);
			set_segment(seg);
		}
	}
}

void
STIRAcquisitionData::xapyb(
	const DataContainer& a_x, const void* ptr_a,
	const DataContainer& a_y, const DataContainer& a_b)
{
	SIRF_DYNAMIC_CAST(const STIRAcquisitionData, x, a_x);
	SIRF_DYNAMIC_CAST(const STIRAcquisitionData, y, a_y);
	SIRF_DYNAMIC_CAST(const STIRAcquisitionData, b, a_b);
	float a = *static_cast<const float*>(ptr_a);
	int n = get_max_segment_num();
	int nx = x.get_max_segment_num();
	int ny = y.get_max_segment_num();
	int nb = b.get_max_segment_num();
	if (n != nx || n != ny || n != nb)
		throw std::runtime_error("binary_op error: operands sizes differ");
	SegmentBySinogram<float>::full_iterator seg_iter;
	SegmentBySinogram<float>::full_iterator sx_iter;
	SegmentBySinogram<float>::full_iterator sy_iter;
	SegmentBySinogram<float>::full_iterator sb_iter;
        TOF_LOOP
	for (int s = 0; s <= n; ++s)
	{
		SegmentBySinogram<float> seg = get_empty_segment_by_sinogram(s TOF_ARG);
		SegmentBySinogram<float> sx = x.get_segment_by_sinogram(s TOF_ARG);
		SegmentBySinogram<float> sy = y.get_segment_by_sinogram(s TOF_ARG);
		SegmentBySinogram<float> sb = b.get_segment_by_sinogram(s TOF_ARG);
		size_t size_seg = seg.size_all();
		size_t size_sx = sx.size_all();
		size_t size_sy = sy.size_all();
		size_t size_sb = sb.size_all();
		if (size_seg != size_sx || size_seg != size_sy || size_seg != size_sb)
			throw std::runtime_error("binary_op error: operands sizes differ");
		for (seg_iter = seg.begin_all(),
			sx_iter = sx.begin_all(), sy_iter = sy.begin_all(), sb_iter = sb.begin_all();
			seg_iter != seg.end_all(); /*empty*/)
			*seg_iter++ = (*sx_iter++) * a + (*sy_iter++) * (*sb_iter++);
		set_segment(seg);
		if (s != 0) {
			seg = get_empty_segment_by_sinogram(-s TOF_ARG);
			sx = x.get_segment_by_sinogram(-s TOF_ARG);
			sy = y.get_segment_by_sinogram(-s TOF_ARG);
			sb = b.get_segment_by_sinogram(-s TOF_ARG);
			for (seg_iter = seg.begin_all(),
				sx_iter = sx.begin_all(), sy_iter = sy.begin_all(), sb_iter = sb.begin_all();
				seg_iter != seg.end_all(); /*empty*/)
				*seg_iter++ = (*sx_iter++) * a + (*sy_iter++) * (*sb_iter++);
			set_segment(seg);
		}
	}
}

std::unique_ptr<STIRAcquisitionData>
STIRAcquisitionDataInFile::get_subset(const std::vector<int>& views) const
{
	auto ptr_ad = new STIRAcquisitionDataInFile(std::move(_data->get_subset(views)));
//	auto ptr_ad = new STIRAcquisitionDataInMemory(std::move(_data->get_subset(views)));
	return std::unique_ptr<STIRAcquisitionData>(ptr_ad);
}

std::unique_ptr<STIRAcquisitionData>
STIRAcquisitionDataInMemory::get_subset(const std::vector<int>& views) const
{
	auto ptr_ad = new STIRAcquisitionDataInMemory(std::move(_data->get_subset(views)));
	return std::unique_ptr<STIRAcquisitionData>(ptr_ad);
}

void
STIRAcquisitionDataInMemory::init()
{
	STIRAcquisitionDataInFile::init();
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
        // STIR format
        parser.add_start_key("Output File Format Parameters");
        // old SIRF format, still accepted
        parser.add_start_key("OutputFileFormat Parameters");
        parser.add_parsing_key("output file format type", &format_sptr);
        parser.add_stop_key("END");
        parser.parse(format_file.c_str());
        if(is_null_ptr(format_sptr))
            throw std::runtime_error("STIRImageData::write: Parsing of output format file (" + format_file + ") failed.\n"
                                     "STIR should have written a warning message why (if you don't see it, redirect warnings to a file using MessageRedirector)\n"
                                     "(see examples/parameter_files/STIR_output_file_format_xxx.par for help).");
    }
    if(is_null_ptr(format_sptr))
        format_sptr = OutputFileFormat<Image3DF>::default_sptr();

    format_sptr->write_to_file(filename, image);
}

void
STIRImageData::sum(void* ptr) const
{
	float* ptr_s = static_cast<float*>(ptr);
        *ptr_s = (float)data().sum();
}

void
STIRImageData::max(void* ptr) const
{
	float* ptr_s = static_cast<float*>(ptr);
        *ptr_s = (float)data().find_max();
}

void
STIRImageData::min(void* ptr) const
{
	float* ptr_s = static_cast<float*>(ptr);
        *ptr_s = (float)data().find_min();
}

void
STIRImageData::dot(const DataContainer& a_x, void* ptr) const
{
	SIRF_DYNAMIC_CAST(const STIRImageData, x, a_x);
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
	float* ptr_s = static_cast<float*>(ptr);
	*ptr_s = (float)s;
}

void
STIRImageData::axpby(
const void* ptr_a, const DataContainer& a_x,
const void* ptr_b, const DataContainer& a_y)
{
	//add deprecation warning
	STIRImageData::xapyb(a_x, ptr_a, a_y, ptr_b);
}

void
STIRImageData::xapyb(
const DataContainer& a_x, const void* ptr_a,
const DataContainer& a_y, const void* ptr_b)
{
	float a = *static_cast<const float*>(ptr_a);
	float b = *static_cast<const float*>(ptr_b);
	SIRF_DYNAMIC_CAST(const STIRImageData, x, a_x);
	SIRF_DYNAMIC_CAST(const STIRImageData, y, a_y);
        data().xapyb(x.data(), a, y.data(), b);
}

void
STIRImageData::xapyb(
const DataContainer& a_x, const void* ptr_a,
const DataContainer& a_y, const DataContainer& a_b)
{
	float a = *static_cast<const float*>(ptr_a);
	SIRF_DYNAMIC_CAST(const STIRImageData, b, a_b);
	SIRF_DYNAMIC_CAST(const STIRImageData, x, a_x);
	SIRF_DYNAMIC_CAST(const STIRImageData, y, a_y);
#if defined(_MSC_VER) && _MSC_VER < 1900
	Image3DF::full_iterator iter;
	Image3DF::const_full_iterator iter_x;
	Image3DF::const_full_iterator iter_y;
	Image3DF::const_full_iterator iter_b;
#else
	typename Array<3, float>::full_iterator iter;
	typename Array<3, float>::const_full_iterator iter_x;
	typename Array<3, float>::const_full_iterator iter_y;
	typename Array<3, float>::const_full_iterator iter_b;
#endif

	if (size() != x.size() || size() != y.size() || size() != b.size())
		throw std::runtime_error("xapyb error: operands sizes differ");

	for (iter = data().begin_all(),
		iter_b = b.data().begin_all(),
		iter_x = x.data().begin_all(), iter_y = y.data().begin_all();
		iter != data().end_all();
		iter++, iter_x++, iter_y++, iter_b++)
		*iter = a * (*iter_x) + (*iter_b) * (*iter_y);
}

void
STIRImageData::xapyb(
	const DataContainer& a_x, const DataContainer& a_a,
	const DataContainer& a_y, const DataContainer& a_b)
{
	SIRF_DYNAMIC_CAST(const STIRImageData, a, a_a);
	SIRF_DYNAMIC_CAST(const STIRImageData, b, a_b);
	SIRF_DYNAMIC_CAST(const STIRImageData, x, a_x);
	SIRF_DYNAMIC_CAST(const STIRImageData, y, a_y);
        data().xapyb(x.data(), a.data(), y.data(), b.data());
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
	return (float)std::sqrt(s);
}

void
STIRImageData::scale(float s)
{
  data() /= s;
}

void
STIRImageData::unary_op(
	const DataContainer& a_x,
	float (*f)(float)
) {
	SIRF_DYNAMIC_CAST(const STIRImageData, x, a_x);
#if defined(_MSC_VER) && _MSC_VER < 1900
	Image3DF::full_iterator iter;
	Image3DF::const_full_iterator iter_x;
#else
	typename Array<3, float>::full_iterator iter;
	typename Array<3, float>::const_full_iterator iter_x;
#endif

	for (iter = data().begin_all(),
		iter_x = x.data().begin_all();
		iter != data().end_all() &&
		iter_x != x.data().end_all();
		iter++, iter_x++)
		*iter = f(*iter_x);
}

void
STIRImageData::semibinary_op(
	const DataContainer& a_x,
	float y, 
	float (*f)(float, float)
){
	SIRF_DYNAMIC_CAST(const STIRImageData, x, a_x);
#if defined(_MSC_VER) && _MSC_VER < 1900
	Image3DF::full_iterator iter;
	Image3DF::const_full_iterator iter_x;
#else
	typename Array<3, float>::full_iterator iter;
	typename Array<3, float>::const_full_iterator iter_x;
#endif

	for (iter = data().begin_all(),
		iter_x = x.data().begin_all();
		iter != data().end_all() &&
		iter_x != x.data().end_all();
		iter++, iter_x++)
		*iter = f(*iter_x, y);
}

void
STIRImageData::binary_op(
	const DataContainer& a_x,
	const DataContainer& a_y,
	float (*f)(float, float)
) {
	SIRF_DYNAMIC_CAST(const STIRImageData, x, a_x);
	SIRF_DYNAMIC_CAST(const STIRImageData, y, a_y);
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
		*iter = f(*iter_x, *iter_y);
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
	SIRF_DYNAMIC_CAST(const Voxels3DF, voxels, *_data);
	CartesianCoordinate3D<float> vs = voxels.get_voxel_size();
	for (int i = 0; i < 3; i++)
		vsize[i] = vs[i + 1];
}

void
STIRImageData::get_data(float* data) const
{
	Coordinate3D<int> min_indices;
	Coordinate3D<int> max_indices;
	if (!_data->get_regular_range(min_indices, max_indices))
		throw LocalisedException("irregular STIR image", __FILE__, __LINE__);
	//std::cout << "trying new const iterator...\n";
	STIRImageData::Iterator_const iter(begin());
	for (int i = 0; iter != end(); ++i, ++iter)
		data[i] = *iter;
	//std::copy(begin(), end(), data);
	//std::copy(image.begin_all(), image.end_all(), data);
	//auto iter = image.begin_all();
	//for (int i = 0; iter != image.end_all(); i++, iter++)
	//	data[i] = *iter;
}

void
STIRImageData::set_data(const float* data)
{
	Image3DF& image = *_data;
	Coordinate3D<int> min_indices;
	Coordinate3D<int> max_indices;
	if (!image.get_regular_range(min_indices, max_indices))
		throw LocalisedException("irregular STIR image", __FILE__, __LINE__);
	size_t n = 1;
	for (int i = 0; i < 3; i++)
		n *= (max_indices[i + 1] - min_indices[i + 1] + 1);
	//std::cout << "trying new iterator...\n";
	STIRImageData::Iterator iter(begin());
	for (int i = 0; iter != end(); ++i, ++iter)
		*iter = data[i];
	//std::copy(data, data + n, begin());
	//std::copy(data, data + n, image.begin_all());
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
    // Need to modify the geom info after changing size
    set_up_geom_info();
}

void
STIRImageData::
zoom_image(const Coord3DF &zooms, const Coord3DF &offsets_in_mm,
           const Coord3DI &new_sizes_in, const stir::ZoomOptions zoom_options)
{
    // We need the underyling image as a VoxelsOnCartesianGrid
    SIRF_DYNAMIC_CAST(Voxels3DF, voxels, this->data());

    int dim[3];
    this->get_dimensions(dim);

    // If any sizes have been set to <= 0, set to image size
    Coord3DI new_sizes(new_sizes_in);
    for (unsigned i=0; i<3; ++i)
        if (new_sizes.at(int(i+1))<=0)
            new_sizes.at(int(i+1)) = dim[i];

    // Zoom the image
    voxels = stir::zoom_image(voxels, zooms, offsets_in_mm, new_sizes, zoom_options);

    // Need to modify the geom info after changing size
    set_up_geom_info();
}

void
STIRImageData::
move_to_scanner_centre(const STIRAcquisitionData &)
{
    this->_data->set_origin(CartesianCoordinate3D<float>{0.f,0.f,0.f});

    // Need to modify the geom info after mod
    set_up_geom_info();
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
        axis_direction /= stir::norm(axis_direction.begin(), axis_direction.end());
        for (int dim = 0; dim < 3; dim++)
            direction[dim][axis] = axis_direction[3 - dim];
    }

    // Initialise the geom info shared pointer
    this->set_geom_info(std::make_shared<VoxelisedGeometricalInfo3D>
                (offset,spacing,size,direction));
}
