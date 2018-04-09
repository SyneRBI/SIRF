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

/*!
\file
\ingroup STIR Extensions
\brief Specification file for data handling types not present in STIR.

\author Evgueni Ovtchinnikov
\author CCP PETMR
*/

#ifndef STIR_DATA_CONTAINER_TYPES
#define STIR_DATA_CONTAINER_TYPES

#include <stdlib.h>

#include <chrono>
#include <fstream>

#include "cstir_shared_ptr.h"
#include "data_handle.h"
#include "stir_types.h"
#include "SIRF/common/data_container.h"

#include "stir/SSRB.h"

using stir::shared_ptr;

class SIRFUtilities {
public:
	static long long milliseconds()
	{
		auto now = std::chrono::system_clock::now();
		auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch());
		return (long long)ms.count();
	}
	static std::string scratch_file_name()
	{
		static int calls = 0;
		char buff[32];
		long long int ms = milliseconds();
		calls++;
		sprintf(buff, "tmp_%d_%lld", calls, ms);
		return std::string(buff);
	}
};

/*!
\ingroup STIR Extensions
\brief STIR ProjDataInterfile wrapper with additional file managing features.

This derived class has additional capability of deleting the file it handles
when an object of this class goes out of existence.
*/

class ProjDataFile : public ProjDataInterfile {
public:
	ProjDataFile(const ProjData& pd, const std::string& filename, bool owns_file = true) :
		ProjDataInterfile(pd.get_exam_info_sptr(),
		pd.get_proj_data_info_sptr(),
		filename, std::ios::in | std::ios::out | std::ios::trunc),
		_filename(filename),
		_owns_file(owns_file)
	{}
	ProjDataFile(shared_ptr<ExamInfo> sptr_exam_info, 
		shared_ptr<ProjDataInfo> sptr_proj_data_info, 
		const std::string& filename, bool owns_file = true) :
		ProjDataInterfile(sptr_exam_info, sptr_proj_data_info,
		filename, std::ios::in | std::ios::out | std::ios::trunc),
		_filename(filename),
		_owns_file(owns_file)
	{}
	~ProjDataFile()
	{
		close_stream();
		clear_stream();
		if (!_owns_file)
			return;
		int err;
		err = std::remove((_filename + ".hs").c_str());
		if (err)
			std::cout << "deleting " << _filename << ".hs "
			<< "failed, please delete manually" << std::endl;
		err = std::remove((_filename + ".s").c_str());
		if (err)
			std::cout << "deleting " << _filename << ".s "
			<< "failed, please delete manually" << std::endl;
	}
	shared_ptr<std::iostream> sino_stream_sptr()
	{
		return sino_stream;
	}
	void close_stream()
	{
		((std::fstream*)sino_stream.get())->close();
	}
	void clear_stream()
	{
		((std::fstream*)sino_stream.get())->clear();
	}
private:
	bool _owns_file;
	std::string _filename;
};

/*!
\ingroup STIR Extensions
\brief STIR ProjData wrapper with added functionality.

This class enjoys some features of STIR ProjData and, additionally,
implements the linear algebra functionality specified by the
abstract base class aDatacontainer, and provides means for the data
storage mode (file/memory) selection.
*/

class PETAcquisitionData : public aDataContainer < float > {
public:
	virtual ~PETAcquisitionData() {}

	// virtual constructors
	virtual PETAcquisitionData* same_acquisition_data
		(shared_ptr<ExamInfo> sptr_exam_info,
		shared_ptr<ProjDataInfo> sptr_proj_data_info) = 0;
	virtual PETAcquisitionData* same_acquisition_data(const ProjData& pd) = 0;
	virtual PETAcquisitionData* same_acquisition_data
		(shared_ptr<ExamInfo> sptr_ei, std::string scanner_name,
		int span = 1, int max_ring_diff = -1, int view_mash_factor = 1) = 0;
	virtual shared_ptr<PETAcquisitionData> new_acquisition_data() = 0;
	virtual aDataContainer<float>* new_data_container() = 0;

	shared_ptr<PETAcquisitionData> single_slice_rebinned_data(
		const int num_segments_to_combine,
		const int num_views_to_combine = 1,
		const int num_tang_poss_to_trim = 0,
		const bool do_normalisation = true,
		const int max_in_segment_num_to_process = -1
		)
	{
		shared_ptr<ProjDataInfo> out_proj_data_info_sptr(
			SSRB(*data()->get_proj_data_info_ptr(),
			num_segments_to_combine,
			num_views_to_combine,
			num_tang_poss_to_trim,
			max_in_segment_num_to_process
			));
		shared_ptr<PETAcquisitionData> sptr(same_acquisition_data
			(data()->get_exam_info_sptr(), out_proj_data_info_sptr));
		SSRB(*sptr, *data(), do_normalisation);
	}

	static std::string storage_scheme()
	{
		static bool initialized = false;
		if (!initialized) {
			_storage_scheme = "file";
			initialized = true;
		}
		return _storage_scheme;
	}

	shared_ptr<ProjData> data()
	{
		return _data;
	}
	const shared_ptr<ProjData> data() const
	{
		return _data;
	}
	void set_data(shared_ptr<ProjData> data)
	{
		_data = data;
	}

	// data import/export
	void fill(float v) { data()->fill(v); }
	void fill(PETAcquisitionData& ad)
	{
		shared_ptr<ProjData> sptr = ad.data();
		data()->fill(*sptr);
	}
	void fill_from(const float* d) { data()->fill_from(d); }
	void copy_to(float* d) { data()->copy_to(d); }
	void write(const char* filename)
	{
		ProjDataFile pd(*data(), filename, false);
		pd.fill(*data());
	}

	// data container methods
	unsigned int items() { return 1; }
	float norm();
	float dot(const aDataContainer<float>& x);
	void mult(float a, const aDataContainer<float>& x);
	void inv(float a, const aDataContainer<float>& x);
	void axpby(float a, const aDataContainer<float>& x,
		float b, const aDataContainer<float>& y);

	// ProjData methods
	int get_num_tangential_poss()
	{
		return data()->get_num_tangential_poss();
	}
	int get_num_views()
	{
		return data()->get_num_views();
	}
	int get_num_sinograms()
	{
		return data()->get_num_sinograms();
	}
	int get_max_segment_num() const
	{
		return data()->get_max_segment_num();
	}
	SegmentBySinogram<float>
		get_segment_by_sinogram(const int segment_num) const
	{
		return data()->get_segment_by_sinogram(segment_num);
	}
	SegmentBySinogram<float>
		get_empty_segment_by_sinogram(const int segment_num) const
	{
		return data()->get_empty_segment_by_sinogram(segment_num);
	}
	virtual Succeeded set_segment(const SegmentBySinogram<float>& s)
	{
		return data()->set_segment(s);
	}
	shared_ptr<ExamInfo> get_exam_info_sptr() const
	{
		return data()->get_exam_info_sptr();
	}
	shared_ptr<ProjDataInfo> get_proj_data_info_sptr() const
	{
		return data()->get_proj_data_info_sptr();
	}

	// ProjData casts
	operator ProjData&() { return *data(); }
	operator const ProjData&() const { return *data(); }
	operator shared_ptr<ProjData>() { return data(); }

protected:
	static std::string _storage_scheme;
	static shared_ptr<PETAcquisitionData> _template;
	shared_ptr<ProjData> _data;

	static shared_ptr<ProjDataInfo> 
		proj_data_info_from_scanner(std::string scanner_name,
		int span = 1, int max_ring_diff = -1, int view_mash_factor = 1)
	{
		shared_ptr<Scanner> sptr_s(Scanner::get_scanner_from_name(scanner_name));
		//std::cout << "scanner: " << sptr_s->get_name().c_str() << '\n';
		if (boost::iequals(sptr_s->get_name(), "unknown")) {
			throw LocalisedException("Unknown scanner", __FILE__, __LINE__);
		}
		int num_views = sptr_s->get_num_detectors_per_ring() / 2 / view_mash_factor;
		int num_tang_pos = sptr_s->get_max_num_non_arccorrected_bins();
		if (max_ring_diff < 0)
			max_ring_diff = sptr_s->get_num_rings() - 1;
		return std::move(ProjDataInfo::construct_proj_data_info
			(sptr_s, span, max_ring_diff, num_views, num_tang_pos, false));
	}
};

/*!
\ingroup STIR Extensions
\brief In-file implementation of PETAcquisitionData.

*/

class PETAcquisitionDataInFile : public PETAcquisitionData {
public:
	PETAcquisitionDataInFile() : _owns_file(false) {}
	PETAcquisitionDataInFile(const char* filename) : _owns_file(false)
	{
		_data = ProjData::read_from_file(filename);
	}
	PETAcquisitionDataInFile(shared_ptr<ExamInfo> sptr_exam_info,
		shared_ptr<ProjDataInfo> sptr_proj_data_info)
	{
		_data.reset(new ProjDataFile
			(sptr_exam_info, sptr_proj_data_info, 
			_filename = SIRFUtilities::scratch_file_name()));
	}
	PETAcquisitionDataInFile(const ProjData& pd) : _owns_file(true)
	{
		_data.reset(new ProjDataFile
		(pd, _filename = SIRFUtilities::scratch_file_name()));
	}
	PETAcquisitionDataInFile
		(shared_ptr<ExamInfo> sptr_ei, std::string scanner_name,
		int span = 1, int max_ring_diff = -1, int view_mash_factor = 1)
	{
		shared_ptr<ProjDataInfo> sptr_pdi =
			PETAcquisitionData::proj_data_info_from_scanner
			(scanner_name, span, max_ring_diff, view_mash_factor);
		ProjDataFile* ptr = new ProjDataFile(sptr_ei, sptr_pdi,
			_filename = SIRFUtilities::scratch_file_name());
		ptr->fill(0.0f);
		_data.reset(ptr);
	}
	shared_ptr<PETAcquisitionData> new_acquisition_data(std::string filename)
	{
		shared_ptr<PETAcquisitionDataInFile> sptr_ad(new PETAcquisitionDataInFile);
		sptr_ad->_data.reset(new ProjDataFile(*data(), filename, false));
		return sptr_ad;
	}

	static void init() {
		static bool initialized = false;
		if (!initialized) {
			_storage_scheme = "file";
			_template.reset(new PETAcquisitionDataInFile());
			initialized = true;
			PETAcquisitionData::storage_scheme();
		}
	}
	static void set_as_template()
	{
		init();
		_storage_scheme = "file";
		_template.reset(new PETAcquisitionDataInFile);
	}

	virtual PETAcquisitionData* same_acquisition_data
		(shared_ptr<ExamInfo> sptr_exam_info,
		shared_ptr<ProjDataInfo> sptr_proj_data_info)
	{
		PETAcquisitionData* ptr_ad = 
			new PETAcquisitionDataInFile(sptr_exam_info, sptr_proj_data_info);
		return ptr_ad;
	}
	virtual PETAcquisitionData* same_acquisition_data(const ProjData& pd)
	{
		PETAcquisitionData* ptr_ad = new PETAcquisitionDataInFile(pd);
		return ptr_ad;
	}
	virtual PETAcquisitionData* same_acquisition_data
		(shared_ptr<ExamInfo> sptr_ei, std::string scanner_name,
		int span = 1, int max_ring_diff = -1, int view_mash_factor = 1)
	{
		PETAcquisitionData* ptr_ad = new PETAcquisitionDataInFile
			(sptr_ei, scanner_name, span, max_ring_diff, view_mash_factor);
		return ptr_ad;
	}
	virtual shared_ptr<PETAcquisitionData> new_acquisition_data()
	{
		init();
		return shared_ptr<PETAcquisitionData>
			(_template->same_acquisition_data(*data()));
	}
	virtual aDataContainer<float>* new_data_container()
	{
		init();
		return (aDataContainer<float>*)_template->same_acquisition_data(*data());
	}

private:
	bool _owns_file;
	std::string _filename;
};

/*!
\ingroup STIR Extensions
\brief In-memory implementation of PETAcquisitionData.

*/

class PETAcquisitionDataInMemory : public PETAcquisitionData {
public:
	PETAcquisitionDataInMemory() {}
	PETAcquisitionDataInMemory(shared_ptr<ExamInfo> sptr_exam_info,
		shared_ptr<ProjDataInfo> sptr_proj_data_info)
	{
		_data = shared_ptr<ProjData>
			(new ProjDataInMemory(sptr_exam_info, sptr_proj_data_info));
	}
	PETAcquisitionDataInMemory(const ProjData& pd)
	{
		_data = shared_ptr<ProjData>
			(new ProjDataInMemory(pd.get_exam_info_sptr(),
			pd.get_proj_data_info_sptr()));
	}
	PETAcquisitionDataInMemory
		(shared_ptr<ExamInfo> sptr_ei, std::string scanner_name,
		int span = 1, int max_ring_diff = -1, int view_mash_factor = 1)
	{
		shared_ptr<ProjDataInfo> sptr_pdi =
			PETAcquisitionData::proj_data_info_from_scanner
			(scanner_name, span, max_ring_diff, view_mash_factor);
		ProjDataInMemory* ptr = new ProjDataInMemory(sptr_ei, sptr_pdi);
		ptr->fill(0.0f);
		_data.reset(ptr);
	}

	static void init() { PETAcquisitionDataInFile::init(); }
	static void set_as_template()
	{
		init();
		_storage_scheme = "memory";
		_template.reset(new PETAcquisitionDataInMemory);
	}

	virtual PETAcquisitionData* same_acquisition_data
		(shared_ptr<ExamInfo> sptr_exam_info,
		shared_ptr<ProjDataInfo> sptr_proj_data_info)
	{
		PETAcquisitionData* ptr_ad =
			new PETAcquisitionDataInMemory(sptr_exam_info, sptr_proj_data_info);
		return ptr_ad;
	}
	virtual PETAcquisitionData* same_acquisition_data(const ProjData& pd)
	{
		PETAcquisitionData* ptr_ad = new PETAcquisitionDataInMemory(pd);
		return ptr_ad;
	}
	virtual PETAcquisitionData* same_acquisition_data
		(shared_ptr<ExamInfo> sptr_ei, std::string scanner_name,
		int span = 1, int max_ring_diff = -1, int view_mash_factor = 1)
	{
		PETAcquisitionData* ptr_ad = new PETAcquisitionDataInMemory
			(sptr_ei, scanner_name, span, max_ring_diff, view_mash_factor);
		return ptr_ad;
	}
	virtual shared_ptr<PETAcquisitionData> new_acquisition_data()
	{
		init();
		return shared_ptr<PETAcquisitionData>
			(_template->same_acquisition_data(*data()));
	}
	virtual aDataContainer<float>* new_data_container()
	{
		init();
		return _template->same_acquisition_data(*data());
	}

};

/*!
\ingroup STIR Extensions
\brief STIR DiscretisedDensity<3, float> wrapper with added functionality.

This class enjoys some features of STIR DiscretisedDensity<3, float> and, 
additioanally, implements the linear algebra functionality specified by the
abstract base class aDatacontainer.
*/

class PETImageData : public aDataContainer<float> {
public:
	PETImageData(){}
	PETImageData(const PETImageData& image)
	{
		_data.reset(image.data().clone());
	}
	PETImageData(const PETAcquisitionData& ad)
	{
		_data.reset(new Voxels3DF(*ad.get_proj_data_info_sptr()));
	}
	PETImageData(const Image3DF& image)
	{
		_data.reset(image.clone());
	}
	PETImageData(const Voxels3DF& v)
	{
		_data.reset(v.clone());
	}
	PETImageData(const ProjDataInfo& pdi)
	{
		_data.reset(new Voxels3DF(pdi));
	}
	PETImageData(shared_ptr<Image3DF> ptr)
	{
		_data = ptr;
	}
	PETImageData(std::string filename)
	{
		_data = read_from_file<Image3DF>(filename);
	}
	PETImageData* same_image_data()
	{
		PETImageData* ptr_image = new PETImageData;
		ptr_image->_data.reset(_data->get_empty_copy());
		return ptr_image;
	}
	shared_ptr<PETImageData> new_image_data()
	{
		return shared_ptr<PETImageData>(same_image_data());
	}
	aDataContainer<float>* new_data_container()
	{
		return same_image_data();
	}
	unsigned int items()
	{
		return 1;
	}
	float norm();
	float dot(const aDataContainer<float>& other);
	void mult(float a, const aDataContainer<float>& x);
	void axpby(float a, const aDataContainer<float>& x,
		float b, const aDataContainer<float>& y);
	Image3DF& data()
	{
		return *_data;
	}
	const Image3DF& data() const
	{
		return *_data;
	}
	Image3DF* data_ptr()
	{
		return _data.get();
	}
	const Image3DF* data_ptr() const
	{
		return _data.get();
	}
	shared_ptr<Image3DF> data_sptr()
	{
		return _data;
	}
	void fill(float v)
	{
		_data->fill(v);
	}
	int get_dimensions(int* dim) const;
	void get_voxel_sizes(float* vsizes) const;
	int get_data(float* data) const;
	int set_data(const float* data);

protected:
	shared_ptr<Image3DF> _data;
};

#endif
