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

#ifndef STIR_DATA_CONTAINER_TYPES
#define STIR_DATA_CONTAINER_TYPES

#include <stdlib.h>

#include <chrono>
#include <fstream>

#include "cstir_shared_ptr.h"
#include "data_handle.h"
#include "stir_types.h"
#include "SIRF/common/data_container.h"

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

//template <typename T>
//class aDataContainer {
//public:
//	virtual ~aDataContainer() {}
//	virtual aDataContainer<T>* new_data_container() = 0;
//	//virtual shared_ptr<aDataContainer<T> > new_data_container() = 0;
//	virtual unsigned int items() = 0;
//	virtual float norm() = 0;
//	virtual T dot(const aDataContainer<T>& dc) = 0;
//	virtual void mult(T a, const aDataContainer<T>& x) = 0;
//	virtual void axpby(
//		T a, const aDataContainer<T>& x,
//		T b, const aDataContainer<T>& y) = 0;
//};

class ProjDataFile : public ProjDataInterfile {
public:
	ProjDataFile(const ProjData& pd, const std::string& filename) :
		ProjDataInterfile(pd.get_exam_info_sptr(),
		pd.get_proj_data_info_sptr(),
		filename, std::ios::in | std::ios::out | std::ios::trunc)
	{}
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
};

class PETAcquisitionData : public aDataContainer < float > {
public:
	virtual ~PETAcquisitionData() {}
	virtual PETAcquisitionData* same_acquisition_data(const ProjData& pd) = 0;
	virtual shared_ptr<PETAcquisitionData> new_acquisition_data() = 0;
	virtual aDataContainer<float>* new_data_container() = 0;
	//virtual shared_ptr<aDataContainer<float> > new_data_container() = 0;
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

	virtual void clear_stream() = 0;
	virtual void close_stream() = 0;

	// ProjData casts
	operator ProjData&() { return *data(); }
	operator const ProjData&() const { return *data(); }

protected:
	//static std::string _storage_scheme;
	static shared_ptr<PETAcquisitionData> _template;
	shared_ptr<ProjData> _data;
};

class PETAcquisitionDataInFile : public PETAcquisitionData {
public:
	PETAcquisitionDataInFile() : _owns_file(false) {}
	PETAcquisitionDataInFile(const char* filename) : _owns_file(false)
	{
		_data = ProjData::read_from_file(filename);
	}
	PETAcquisitionDataInFile(const ProjData& pd) : _owns_file(true)
	{
		_data.reset(new ProjDataFile
		(pd, _filename = SIRFUtilities::scratch_file_name()));
	}
	~PETAcquisitionDataInFile()
	{
		_data.reset();
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

	static void init() {
		static bool initialized = false;
		if (!initialized) {
			_template.reset(new PETAcquisitionDataInFile());
			initialized = true;
		}
	}
	static void set_as_template()
	{
		init();
		_template.reset(new PETAcquisitionDataInFile);
	}

	PETAcquisitionData* same_acquisition_data(const ProjData& pd)
	{
		PETAcquisitionData* ptr_ad = new PETAcquisitionDataInFile(pd);
		return ptr_ad;
	}
	shared_ptr<PETAcquisitionData> new_acquisition_data()
	{
		init();
		return shared_ptr<PETAcquisitionData>
			(_template->same_acquisition_data(*data()));
	}
	//shared_ptr<aDataContainer<float> > new_data_container()
	aDataContainer<float>* new_data_container()
	{
		init();
		return (aDataContainer<float>*)_template->same_acquisition_data(*data());
		//return shared_ptr<aDataContainer<float> >
		//	(_template->same_acquisition_data(*data()));
	}

	void clear_stream()
	{
		((ProjDataFile*)_data.get())->clear_stream();
	}
	void close_stream()
	{
		((ProjDataFile*)_data.get())->close_stream();
	}

private:
	bool _owns_file;
	std::string _filename;
};

class PETAcquisitionDataInMemory : public PETAcquisitionData {
public:
	PETAcquisitionDataInMemory() {}
	PETAcquisitionDataInMemory(const ProjData& pd)
	{
		_data = shared_ptr<ProjData>
			(new ProjDataInMemory(pd.get_exam_info_sptr(),
			pd.get_proj_data_info_sptr()));
	}

	static void init() { PETAcquisitionDataInFile::init(); }
	static void set_as_template()
	{
		init();
		_template.reset(new PETAcquisitionDataInMemory);
	}

	PETAcquisitionData* same_acquisition_data(const ProjData& pd)
	{
		PETAcquisitionData* ptr_ad = new PETAcquisitionDataInMemory(pd);
		return ptr_ad;
	}
	shared_ptr<PETAcquisitionData> new_acquisition_data()
	{
		init();
		return shared_ptr<PETAcquisitionData>
			(_template->same_acquisition_data(*data()));
	}
	aDataContainer<float>* new_data_container()
	{
		init();
		return _template->same_acquisition_data(*data());
	}
	//shared_ptr<aDataContainer<float> > new_data_container()
	//{
	//	init();
	//	return shared_ptr<aDataContainer<float> >
	//		(_template->same_acquisition_data(*data()));
	//}

	void clear_stream() {}
	void close_stream()	{}
};

class PETImageData : public aDataContainer<float> {
public:
	PETImageData(){}
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
	//shared_ptr<aDataContainer<float> > new_data_container()
	//{
	//	return shared_ptr<aDataContainer<float> >(same_image_data());
	//}
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

protected:
	shared_ptr<Image3DF> _data;
};

#endif
