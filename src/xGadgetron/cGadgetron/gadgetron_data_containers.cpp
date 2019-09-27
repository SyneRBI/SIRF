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

/*!
\file
\ingroup Gadgetron Data Containers
\brief Implementation file for SIRF data container classes for Gadgetron data.

\author Evgueni Ovtchinnikov
\author CCP PETMR
*/
#include <cmath>
#include <iomanip>

#include "sirf/Gadgetron/cgadgetron_shared_ptr.h"
#include "sirf/Gadgetron/gadgetron_data_containers.h"

using namespace gadgetron;
using namespace sirf;

std::string MRAcquisitionData::_storage_scheme;
shared_ptr<MRAcquisitionData> MRAcquisitionData::acqs_templ_;

static std::string get_date_time_string()
{
    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );

    std::stringstream str;
    str << timeinfo->tm_year+1900 << "-"
        << std::setw(2) << std::setfill('0') << timeinfo->tm_mon+1 << "-"
        << std::setw(2) << std::setfill('0') << timeinfo->tm_mday << " "
        << std::setw(2) << std::setfill('0') << timeinfo->tm_hour << ":"
        << std::setw(2) << std::setfill('0') << timeinfo->tm_min << ":"
        << std::setw(2) << std::setfill('0') << timeinfo->tm_sec;
    return str.str();
}

void 
MRAcquisitionData::write(const std::string &filename) const
{
	Mutex mtx;
	mtx.lock();
	shared_ptr<ISMRMRD::Dataset> dataset
		(new ISMRMRD::Dataset(filename.c_str(), "/dataset", true));
	dataset->writeHeader(acqs_info_);
	mtx.unlock();
	int n = number();
	ISMRMRD::Acquisition a;
	for (int i = 0; i < n; i++) {
		get_acquisition(i, a);
		//std::cout << i << ' ' << a.idx().repetition << '\n';
		//if (TO_BE_IGNORED(a)) {
		//	continue;
		//}
		mtx.lock();
		dataset->appendAcquisition(a);
		mtx.unlock();
	}
}

void
MRAcquisitionData::read( const std::string& filename_ismrmrd_with_ext )
{
	
	bool const verbose = true;

	if( verbose )
		std::cout<< "Started reading acquisitions from " << filename_ismrmrd_with_ext << std::endl;
	try
	{

		ISMRMRD::Dataset d(filename_ismrmrd_with_ext.c_str(),"dataset", false);

		d.readHeader(this->acqs_info_);

		uint32_t num_acquis = d.getNumberOfAcquisitions();
		for( uint32_t i_acqu=0; i_acqu<num_acquis; i_acqu++)
		{
			if( verbose )
			{
				if( i_acqu%( num_acquis/10 ) == 0 )
					std::cout << std::ceil( float(i_acqu)/num_acquis*100 )<< " % " << " done."<< std::endl;
			}

			ISMRMRD::Acquisition acq;
			d.readAcquisition( i_acqu, acq);

			if( TO_BE_IGNORED(acq) )
				continue;
			else
				this->append_acquisition( acq );
		}
		if( verbose )
			std::cout<< "Finished reading acquisitions from " << filename_ismrmrd_with_ext << std::endl;
	}
	catch( std::runtime_error& e)
	{
		std::cerr << "An exception was caught reading " << filename_ismrmrd_with_ext << std::endl;
		std::cerr << e.what() <<std::endl;
		throw;
	}
}

bool
MRAcquisitionData::undersampled() const
{
	ISMRMRD::IsmrmrdHeader header;
	ISMRMRD::deserialize(acqs_info_.c_str(), header);
	ISMRMRD::Encoding e = header.encoding[0];
	return e.parallelImaging.is_present() &&
		e.parallelImaging().accelerationFactor.kspace_encoding_step_1 > 1;
}

int 
MRAcquisitionData::get_acquisitions_dimensions(size_t ptr_dim) const
{
	ISMRMRD::Acquisition acq;
	int* dim = (int*)ptr_dim;

	int na = number();
	int ms, ns;
	int mc, nc;
	int my, ny;
	int slice = 0;
	int y = 0;
	// assume all dimensions (samples, coils [, acqs per slice]) regular
	int nrd = sorted() ? 3 : 2;
	// number of regular readouts
	int nrr = 0;
	//int not_reg = 0;
	for (; y < na;) {
		for (; y < na && sorted();) {
			get_acquisition(y, acq);
			if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_SLICE))
				break;
			y++;
		}
		if (y >= na)
			break;
		ny = 0;
		for (; y < na; y++) {
			get_acquisition(y, acq);
			if (TO_BE_IGNORED(acq)) // not a regular acquisition
				continue;
			ns = acq.number_of_samples();
			nc = acq.active_channels();
			nrr += ns*nc;
			if (slice == 0) {
				ms = ns;
				mc = nc;
			}
			else {
				if (ms != ns)
					nrd = 0;
				else if (mc != nc && nrd > 1)
					nrd = 1;
			}
			ny++;
			if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE) && sorted())
				break;
		}
		if (slice == 0) {
			my = ny;
		}
		else {
			if (my != ny && nrd > 2) {
				nrd = 2;
			}
			//if (my != ny)
			//	not_reg = 1;
		}
		slice++;
	}

	int reg_size = 1;
	if (nrd > 0) {
		dim[0] = ms;
		reg_size *= ms;
	}
	if (nrd > 1) {
		dim[1] = mc;
		reg_size *= mc;
	}
	if (nrd > 2) {
		dim[2] = my;
		reg_size *= my;
	}
	dim[nrd] = nrr / reg_size;
	return nrd;
}

void
MRAcquisitionData::get_data(complex_float_t* z, int all)
{
	ISMRMRD::Acquisition acq;
	unsigned int na = number();
	unsigned int n = 0;
	for (unsigned int a = 0, i = 0; a < na; a++) {
		get_acquisition(a, acq);
		if (!all && TO_BE_IGNORED(acq)) {
			std::cout << "ignoring acquisition " << a << '\n';
			continue;
		}
		n++;
		unsigned int nc = acq.active_channels();
		unsigned int ns = acq.number_of_samples();
		for (unsigned int c = 0; c < nc; c++) {
			for (unsigned int s = 0; s < ns; s++, i++) {
				z[i] = acq.data(s, c);
			}
		}
	}
}

void 
MRAcquisitionData::axpby
(complex_float_t a, const ISMRMRD::Acquisition& acq_x,
	complex_float_t b, ISMRMRD::Acquisition& acq_y)
{
	complex_float_t* px;
	complex_float_t* py;
	for (px = acq_x.data_begin(), py = acq_y.data_begin();
		px != acq_x.data_end() && py != acq_y.data_end(); px++, py++) {
		if (b == complex_float_t(0.0))
			*py = a*complex_float_t(*px);
		else
			*py = a*complex_float_t(*px) + b*complex_float_t(*py);
	}
}

void
MRAcquisitionData::multiply
(const ISMRMRD::Acquisition& acq_x, ISMRMRD::Acquisition& acq_y)
{
	complex_float_t* px;
	complex_float_t* py;
	for (px = acq_x.data_begin(), py = acq_y.data_begin();
		px != acq_x.data_end() && py != acq_y.data_end(); px++, py++) {
		*py = complex_float_t(*px) * complex_float_t(*py);
	}
}

void
MRAcquisitionData::divide
(const ISMRMRD::Acquisition& acq_x, ISMRMRD::Acquisition& acq_y)
{
	complex_float_t* px;
	complex_float_t* py;
	for (px = acq_x.data_begin(), py = acq_y.data_begin();
		px != acq_x.data_end() && py != acq_y.data_end(); px++, py++) {
		// TODO: check for zero denominator
		*py = complex_float_t(*px) / complex_float_t(*py);
	}
}

complex_float_t
MRAcquisitionData::dot
(const ISMRMRD::Acquisition& acq_a, const ISMRMRD::Acquisition& acq_b)
{
	complex_float_t* pa;
	complex_float_t* pb;
	complex_float_t z = 0;
	for (pa = acq_a.data_begin(), pb = acq_b.data_begin();
		pa != acq_a.data_end() && pb != acq_b.data_end(); pa++, pb++) {
		z += std::conj(*pb) * (*pa);
	}
	return z;
}

float 
MRAcquisitionData::norm(const ISMRMRD::Acquisition& acq_a)
{
	complex_float_t* pa;
	float r = 0;
	for (pa = acq_a.data_begin(); pa != acq_a.data_end(); pa++) {
		complex_float_t z = std::conj(*pa) * (*pa);
		r += z.real();
	}
	r = sqrt(r);
	return r;
}

void
MRAcquisitionData::dot(const DataContainer& dc, void* ptr) const
{
	//MRAcquisitionData& other = (MRAcquisitionData&)dc;
	DYNAMIC_CAST(const MRAcquisitionData, other, dc);
	int n = number();
	int m = other.number();
	complex_float_t z = 0;
	ISMRMRD::Acquisition a;
	ISMRMRD::Acquisition b;
	for (int i = 0, j = 0; i < n && j < m;) {
		get_acquisition(i, a);
		if (TO_BE_IGNORED(a)) {
			i++;
			continue;
		}
		other.get_acquisition(j, b);
		if (TO_BE_IGNORED(b)) {
			j++;
			continue;
		}
		z += MRAcquisitionData::dot(a, b);
		i++;
		j++;
	}
	complex_float_t* ptr_z = (complex_float_t*)ptr;
	*ptr_z = z;
}

void
MRAcquisitionData::axpby(
const void* ptr_a, const DataContainer& a_x,
const void* ptr_b, const DataContainer& a_y)
{
	complex_float_t a = *(complex_float_t*)ptr_a;
	complex_float_t b = *(complex_float_t*)ptr_b;
	DYNAMIC_CAST(const MRAcquisitionData, x, a_x);
	DYNAMIC_CAST(const MRAcquisitionData, y, a_y);
	//MRAcquisitionData& x = (MRAcquisitionData&)a_x;
	//MRAcquisitionData& y = (MRAcquisitionData&)a_y;
	int m = x.number();
	int n = y.number();
	ISMRMRD::Acquisition ax;
	ISMRMRD::Acquisition ay;
	for (int i = 0, j = 0; i < n && j < m;) {
		y.get_acquisition(i, ay);
		x.get_acquisition(j, ax);
		if (TO_BE_IGNORED(ay)) {
			std::cout << i << " ignored (ay)\n";
			i++;
			continue;
		}
		if (TO_BE_IGNORED(ax)) {
			std::cout << j << " ignored (ax)\n";
			j++;
			continue;
		}
		MRAcquisitionData::axpby(a, ax, b, ay);
		append_acquisition(ay);
		i++;
		j++;
	}
}

void
MRAcquisitionData::multiply(
const DataContainer& a_x,
const DataContainer& a_y)
{
	//MRAcquisitionData& x = (MRAcquisitionData&)a_x;
	//MRAcquisitionData& y = (MRAcquisitionData&)a_y;
	DYNAMIC_CAST(const MRAcquisitionData, x, a_x);
	DYNAMIC_CAST(const MRAcquisitionData, y, a_y);
	int m = x.number();
	int n = y.number();
	ISMRMRD::Acquisition ax;
	ISMRMRD::Acquisition ay;
	for (int i = 0, j = 0; i < n && j < m;) {
		y.get_acquisition(i, ay);
		x.get_acquisition(j, ax);
		if (TO_BE_IGNORED(ay)) {
			std::cout << i << " ignored (ay)\n";
			i++;
			continue;
		}
		if (TO_BE_IGNORED(ax)) {
			std::cout << j << " ignored (ax)\n";
			j++;
			continue;
		}
		MRAcquisitionData::multiply(ax, ay);
		append_acquisition(ay);
		i++;
		j++;
	}
}

void
MRAcquisitionData::divide(
const DataContainer& a_x,
const DataContainer& a_y)
{
	//MRAcquisitionData& x = (MRAcquisitionData&)a_x;
	//MRAcquisitionData& y = (MRAcquisitionData&)a_y;
	DYNAMIC_CAST(const MRAcquisitionData, x, a_x);
	DYNAMIC_CAST(const MRAcquisitionData, y, a_y);
	int m = x.number();
	int n = y.number();
	ISMRMRD::Acquisition ax;
	ISMRMRD::Acquisition ay;
	for (int i = 0, j = 0; i < n && j < m;) {
		y.get_acquisition(i, ay);
		x.get_acquisition(j, ax);
		if (TO_BE_IGNORED(ay)) {
			std::cout << i << " ignored (ay)\n";
			i++;
			continue;
		}
		if (TO_BE_IGNORED(ax)) {
			std::cout << j << " ignored (ax)\n";
			j++;
			continue;
		}
		MRAcquisitionData::divide(ax, ay);
		append_acquisition(ay);
		i++;
		j++;
	}
}

float 
MRAcquisitionData::norm() const
{
	int n = number();
	float r = 0;
	ISMRMRD::Acquisition a;
	for (int i = 0; i < n; i++) {
		get_acquisition(i, a);
		if (TO_BE_IGNORED(a)) {
			continue;
		}
		float s = MRAcquisitionData::norm(a);
		r += s*s;
	}
	return sqrt(r);
}

MRAcquisitionData*
MRAcquisitionData::clone_base() const
{
	MRAcquisitionData* ptr_ad =
		acqs_templ_->same_acquisitions_container(this->acqs_info_);
	for (int i = 0; i < number(); i++) {
		ISMRMRD::Acquisition acq;
		get_acquisition(i, acq);
		ptr_ad->append_acquisition(acq);
	}
	return ptr_ad;
}

void
MRAcquisitionData::sort()
{
	typedef std::array<int, 4> tuple;
	int na = number();
	int last = -1;
	int max_phase = 0;
	tuple t;
	std::vector<tuple> vt;
	for (int i = 0; i < na; i++) {
		ISMRMRD::Acquisition acq;
		get_acquisition(i, acq);
		if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_MEASUREMENT))
			last = i;
		t[0] = acq.idx().repetition;
		t[1] = acq.idx().phase;
		t[2] = acq.idx().slice;
		t[3] = acq.idx().kspace_encode_step_1;
		vt.push_back(t);
		if (t[1] > max_phase)
			max_phase = t[1];
	}
	if (last > -1)
		vt[last][1] = max_phase;

	index_.resize(na);

	if( na <= 0 )
		std::cerr << "WARNING: You try to sort an empty container of acquisition data." << std::endl;
	else
		Multisort::sort( vt, &index_[0] );

	sorted_ = true;
}

void
MRAcquisitionData::sort_by_time()
{
	typedef std::array<uint32_t , 1>  tuple;
	tuple t;
	std::vector< tuple > vt;
	size_t const num_acquis = this->number();

	for(size_t i=0; i<num_acquis; i++)
	{
		ISMRMRD::Acquisition acq;
		get_acquisition(i, acq);
		t[0] = acq.acquisition_time_stamp();
		vt.push_back( t );
	}

	index_.resize(num_acquis);
	
	if( num_acquis == 0 )
		std::cerr << "WARNING: You try to sort by time an empty container of acquisition data." << std::endl;
	else
		Multisort::sort( vt, &index_[0] );

}

AcquisitionsFile::AcquisitionsFile
(std::string filename, bool create_file, AcquisitionsInfo info)
{
	own_file_ = create_file;
	filename_ = filename;

	Mutex mtx;
	mtx.lock();
	dataset_ = shared_ptr<ISMRMRD::Dataset>
		(new ISMRMRD::Dataset(filename.c_str(), "/dataset", create_file));
	if (!create_file) {
		dataset_->readHeader(acqs_info_);
	}
	else {
		acqs_info_ = info;
		dataset_->writeHeader(acqs_info_);
	}
	mtx.unlock();
}

AcquisitionsFile::AcquisitionsFile(AcquisitionsInfo info)
{
	own_file_ = true;
	filename_ = xGadgetronUtilities::scratch_file_name();
	Mutex mtx;
	mtx.lock();
	dataset_ = shared_ptr<ISMRMRD::Dataset>
		(new ISMRMRD::Dataset(filename_.c_str(), "/dataset", true));
	acqs_info_ = info;
	dataset_->writeHeader(acqs_info_);
	mtx.unlock();
}

AcquisitionsFile::~AcquisitionsFile() 
{
	dataset_.reset();
	if (own_file_) {
		Mutex mtx;
		mtx.lock();
		std::remove(filename_.c_str());
		mtx.unlock();
	}
}

void 
AcquisitionsFile::take_over(AcquisitionsFile& af)
{
	acqs_info_ = af.acquisitions_info();
	sorted_ = af.sorted();
	index_ = af.index();
	dataset_ = af.dataset_;
	if (own_file_) {
		Mutex mtx;
		mtx.lock();
		std::remove(filename_.c_str());
		mtx.unlock();
	}
	filename_ = af.filename_;
	own_file_ = af.own_file_;
	af.own_file_ = false;
}

unsigned int 
AcquisitionsFile::items() const
{
	Mutex mtx;
	mtx.lock();
	unsigned int na = dataset_->getNumberOfAcquisitions();
	mtx.unlock();
	return na;
}

void 
AcquisitionsFile::get_acquisition(unsigned int num, ISMRMRD::Acquisition& acq) const
{
	int ind = index(num);
	Mutex mtx;
	mtx.lock();
	dataset_->readAcquisition(ind, acq);
	//dataset_->readAcquisition(index(num), acq); // ??? does not work!
	mtx.unlock();
}

void 
AcquisitionsFile::append_acquisition(ISMRMRD::Acquisition& acq)
{
	Mutex mtx;
	mtx.lock();
	dataset_->appendAcquisition(acq);
	mtx.unlock();
}

void 
AcquisitionsFile::copy_acquisitions_info(const MRAcquisitionData& ac)
{
	acqs_info_ = ac.acquisitions_info();
	Mutex mtx;
	mtx.lock();
	dataset_->writeHeader(acqs_info_);
	mtx.unlock();
}

void 
AcquisitionsFile::write_acquisitions_info()
{
	Mutex mtx;
	mtx.lock();
	dataset_->writeHeader(acqs_info_);
	mtx.unlock();
}

void
AcquisitionsFile::set_data(const complex_float_t* z, int all)
{
	AcquisitionsFile ac(acqs_info_);
	ISMRMRD::Acquisition acq;
	int na = number();
	for (int a = 0, i = 0; a < na; a++) {
		get_acquisition(a, acq);
		if (!all && TO_BE_IGNORED(acq)) {
			std::cout << "ignoring acquisition " << a << '\n';
			continue;
		}
		unsigned int nc = acq.active_channels();
		unsigned int ns = acq.number_of_samples();
		for (int c = 0; c < nc; c++)
			for (int s = 0; s < ns; s++, i++)
				acq.data(s, c) = z[i];
		ac.append_acquisition(acq);
	}
	take_over(ac);
}

void
AcquisitionsFile::copy_acquisitions_data(const MRAcquisitionData& ac)
{
	AcquisitionsFile af(acqs_info_);
	ISMRMRD::Acquisition acq;
	int na = number();
	assert(na == ac.number());
	for (int a = 0, i = 0; a < na; a++) {
		ac.get_acquisition(a, acq);
		af.append_acquisition(acq);
	}
	take_over(af);
}

void
AcquisitionsVector::set_data(const complex_float_t* z, int all)
{
	int na = number();
	for (int a = 0, i = 0; a < na; a++) {
		ISMRMRD::Acquisition& acq = *acqs_[a];
		if (!all && TO_BE_IGNORED(acq)) {
			std::cout << "ignoring acquisition " << a << '\n';
			continue;
		}
		unsigned int nc = acq.active_channels();
		unsigned int ns = acq.number_of_samples();
		for (int c = 0; c < nc; c++)
			for (int s = 0; s < ns; s++, i++)
				acq.data(s, c) = z[i];
	}
}

void
AcquisitionsVector::copy_acquisitions_data(const MRAcquisitionData& ac)
{
	ISMRMRD::Acquisition acq_dst;
	ISMRMRD::Acquisition acq_src;
	int na = number();
	assert(na == ac.number());
	for (int a = 0, i = 0; a < na; a++) {
		ac.get_acquisition(a, acq_src);
		ISMRMRD::Acquisition& acq_dst = *acqs_[a];
		unsigned int nc = acq_dst.active_channels();
		unsigned int ns = acq_dst.number_of_samples();
		assert(nc == acq_src.active_channels());
		assert(ns == acq_src.number_of_samples());
		for (int c = 0; c < nc; c++)
			for (int s = 0; s < ns; s++, i++)
				acq_dst.data(s, c) = acq_src.data(s, c);
	}
}

void
GadgetronImageData::dot(const DataContainer& dc, void* ptr) const
{
	//GadgetronImageData& ic = (GadgetronImageData&)dc;
	DYNAMIC_CAST(const GadgetronImageData, ic, dc);
	complex_float_t z = 0;
	for (unsigned int i = 0; i < number() && i < ic.number(); i++) {
		const ImageWrap& u = image_wrap(i);
		const ImageWrap& v = ic.image_wrap(i);
		z += u.dot(v);
	}
	complex_float_t* ptr_z = (complex_float_t*)ptr;
	*ptr_z = z;
}

void
GadgetronImageData::axpby(
const void* ptr_a, const DataContainer& a_x,
const void* ptr_b, const DataContainer& a_y)
{
	complex_float_t a = *(complex_float_t*)ptr_a;
	complex_float_t b = *(complex_float_t*)ptr_b;
	DYNAMIC_CAST(const GadgetronImageData, x, a_x);
	DYNAMIC_CAST(const GadgetronImageData, y, a_y);
	//GadgetronImageData& x = (GadgetronImageData&)a_x;
	//GadgetronImageData& y = (GadgetronImageData&)a_y;
	ImageWrap w(x.image_wrap(0));
	complex_float_t zero(0.0, 0.0);
	complex_float_t one(1.0, 0.0);
	for (unsigned int i = 0; i < x.number() && i < y.number(); i++) {
		const ImageWrap& u = x.image_wrap(i);
		const ImageWrap& v = y.image_wrap(i);
		w.axpby(a, u, zero);
		w.axpby(b, v, one);
		append(w);
	}
}

void
GadgetronImageData::multiply(
const DataContainer& a_x,
const DataContainer& a_y)
{
	//GadgetronImageData& x = (GadgetronImageData&)a_x;
	//GadgetronImageData& y = (GadgetronImageData&)a_y;
	DYNAMIC_CAST(const GadgetronImageData, x, a_x);
	DYNAMIC_CAST(const GadgetronImageData, y, a_y);
	for (unsigned int i = 0; i < x.number() && i < y.number(); i++) {
		ImageWrap w(x.image_wrap(i));
		w.multiply(y.image_wrap(i));
		append(w);
	}
}

void
GadgetronImageData::divide(
const DataContainer& a_x,
const DataContainer& a_y)
{
	//GadgetronImageData& x = (GadgetronImageData&)a_x;
	//GadgetronImageData& y = (GadgetronImageData&)a_y;
	DYNAMIC_CAST(const GadgetronImageData, x, a_x);
	DYNAMIC_CAST(const GadgetronImageData, y, a_y);
	for (unsigned int i = 0; i < x.number() && i < y.number(); i++) {
		ImageWrap w(x.image_wrap(i));
		w.divide(y.image_wrap(i));
		append(w);
	}
}

float 
GadgetronImageData::norm() const
{
	float r = 0;
	for (unsigned int i = 0; i < number(); i++) {
		const ImageWrap& u = image_wrap(i);
		float s = u.norm();
		r += s*s;
	}
	r = sqrt(r);
	return r;
}

void
GadgetronImagesVector::sort()
{
	typedef std::array<float, 3> tuple;
	int ni = number();
	tuple t;
	std::vector<tuple> vt;
	for (int i = 0; i < ni; i++) {
      ImageWrap& iw = image_wrap(i);
      ISMRMRD::ImageHeader& head = iw.head();
		t[0] = head.contrast;
        t[1] = head.repetition;
        // Calculate the projection of the position in the slice direction
        t[2] = head.position[0] * head.slice_dir[0] +
               head.position[1] * head.slice_dir[1] +
               head.position[2] * head.slice_dir[2];
		vt.push_back(t);
#ifndef NDEBUG
        std::cout << "Before sorting. Image " << i << "/" << ni <<  ", Contrast: " << t[0] << ", Repetition: " << t[1] << ", Projection: " << t[2] << "\n";
#endif
	}

	index_.resize(ni);
	Multisort::sort(vt, &index_[0] );
	sorted_ = true;

	// quick fix for the problem of compatibility with image data iterators
	std::vector<gadgetron::shared_ptr<ImageWrap> > sorted_images;
	for (int i = 0; i < ni; i++)
		sorted_images.push_back(sptr_image_wrap(i));
	images_ = sorted_images;
	index_.resize(0);

#ifndef NDEBUG
    std::cout << "After sorting...\n";
    for (int i = 0; i < ni; i++) {
      ImageWrap& iw = image_wrap(i);
      ISMRMRD::ImageHeader& head = iw.head();
		t[0] = head.contrast;
        t[1] = head.repetition;
        // Calculate the projection of the position in the slice direction
        t[2] = head.position[0] * head.slice_dir[0] +
               head.position[1] * head.slice_dir[1] +
               head.position[2] * head.slice_dir[2];
        std::cout << "Image " << i << "/" << ni <<  ", Contrast: " << t[0] << ", Repetition: " << t[1] << ", Projection: " << t[2] << "\n";
	}
#endif
}

std::shared_ptr<std::vector<std::string> >
group_names_sptr(const char* filename)
{
	hid_t    file;
	hid_t    root;
	hid_t    group;
	herr_t   err;
	herr_t   status;
	hsize_t  nobj;
	const int MAX_NAME = 1024;
	char group_name[MAX_NAME];
	char var_name[MAX_NAME];
	std::shared_ptr<std::vector<std::string> >
		sptr_names(new std::vector<std::string>);
	std::vector<std::string>& names = *sptr_names;

	file = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
	if (file < 0)
		return sptr_names;
	root = H5Gopen(file, "/", H5P_DEFAULT);
	if (root < 0) {
		status = H5Fclose(file);
		return sptr_names;
	}
	H5Gget_objname_by_idx(root, 0, group_name, (size_t)MAX_NAME);
	names.push_back(std::string(group_name));
	group = H5Gopen(root, group_name, H5P_DEFAULT);
	err = H5Gget_num_objs(group, &nobj);
	//printf("%d objects\n", nobj);
	for (int i = 0; i < nobj; i++) {
		H5Gget_objname_by_idx(group, (hsize_t)i,
			var_name, (size_t)MAX_NAME);
		names.push_back(std::string(var_name));
	}

	H5Gclose(group);
	H5Gclose(root);
	status = H5Fclose(file);

	return sptr_names;
}

int
GadgetronImageData::read(std::string filename) 
{
	std::shared_ptr<std::vector<std::string> > sptr_names;
	sptr_names = group_names_sptr(filename.c_str());
	std::vector<std::string>& names = *sptr_names;
	int ng = names.size();
	const char* group = names[0].c_str();
	printf("group %s\n", group);
	for (int i = 0; i < ng; i++) {
		const char* var = names[i].c_str();
		if (!i)
			continue;

		printf("variable %s\n", var);
		ISMRMRD::ISMRMRD_Dataset dataset;
		ISMRMRD::ISMRMRD_Image im;
		ismrmrd_init_dataset(&dataset, filename.c_str(), group);
		ismrmrd_open_dataset(&dataset, false);
        // Read XML meta data. Can't use readHeader as this throws error if no XML present
        this->acqs_info_ = ismrmrd_read_header(&dataset);
		int num_im = ismrmrd_get_number_of_images(&dataset, var);
		std::cout << "number of images: " << num_im << '\n';
		ismrmrd_init_image(&im);
		ismrmrd_read_image(&dataset, var, 0, &im);
		printf("image data type: %d\n", im.head.data_type);
		ismrmrd_cleanup_image(&im);
		ismrmrd_close_dataset(&dataset);


		shared_ptr<ISMRMRD::Dataset> sptr_dataset
			(new ISMRMRD::Dataset(filename.c_str(), group, false));

		for (int i = 0; i < num_im; i++) {
			shared_ptr<ImageWrap> sptr_iw(new ImageWrap(im.head.data_type, *sptr_dataset, var, i));
			//sptr_iw->read(*sptr_dataset, var, i);
			append(*sptr_iw);
			//images_.push_back(sptr_iw);
		}
		//int dim[3];
		//sptr_iw->get_dim(dim);
		//std::cout << "image dimensions: "
		//	<< dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';
	}

    this->set_up_geom_info();
	return 0;
}

void
GadgetronImageData::write(const std::string &filename, const std::string &groupname) const
{
	//if (images_.size() < 1)
	if (number() < 1)
		return;
    // If the groupname hasn't been set, use the current date and time.
    std::string group = groupname;
    if (group.empty())
        group = get_date_time_string();
	Mutex mtx;
	mtx.lock();
	ISMRMRD::Dataset dataset(filename.c_str(), group.c_str());
    dataset.writeHeader(acqs_info_);
	mtx.unlock();
	for (unsigned int i = 0; i < number(); i++) {
		const ImageWrap& iw = image_wrap(i);
		iw.write(dataset);
	}
}

void
GadgetronImageData::get_data(complex_float_t* data) const
{
	int dim[4];
	for (unsigned int i = 0; i < number(); i++) {
		const ImageWrap& iw = image_wrap(i);
		size_t n = iw.get_dim(dim);
		iw.get_complex_data(data);
		data += n;
	}
}

void
GadgetronImageData::set_data(const complex_float_t* z)
{
	int dim[4];
	for (unsigned int i = 0; i < number(); i++) {
		ImageWrap& iw = image_wrap(i);
		size_t n = iw.get_dim(dim);
		iw.set_complex_data(z);
		z += n;
	}
}

void
GadgetronImageData::get_real_data(float* data) const
{
	int dim[4];
	for (unsigned int i = 0; i < number(); i++) {
		const ImageWrap& iw = image_wrap(i);
		size_t n = iw.get_dim(dim);
		iw.get_data(data);
		data += n;
	}
}

void
GadgetronImageData::set_real_data(const float* z)
{
	int dim[4];
	for (unsigned int i = 0; i < number(); i++) {
		ImageWrap& iw = image_wrap(i);
		size_t n = iw.get_dim(dim);
		iw.set_data(z);
		z += n;
	}
}

GadgetronImagesVector::GadgetronImagesVector
(const GadgetronImagesVector& images) :
images_()
{
	DYNAMIC_CAST(const GadgetronImageData, imgs, images);
	set_meta_data(imgs.get_meta_data());
	for (unsigned int i = 0; i < images.number(); i++) {
		const ImageWrap& u = images.image_wrap(i);
		append(u);
	}
    this->set_up_geom_info();
}

GadgetronImagesVector::GadgetronImagesVector
(GadgetronImagesVector& images, const char* attr, const char* target) : 
images_()
{
	DYNAMIC_CAST(const GadgetronImageData, imgs, images);
	set_meta_data(imgs.get_meta_data());
	for (unsigned int i = 0; i < images.number(); i++) {
		const ImageWrap& u = images.image_wrap(i);
		std::string atts = u.attributes();
		ISMRMRD::MetaContainer mc;
		ISMRMRD::deserialize(atts.c_str(), mc);
		size_t l = mc.length(attr);
		std::string value;
		for (int j = 0; j < l; j++) {
			if (j)
				value += " ";
			value += mc.as_str(attr, j);
		}
		//std::cout << value.c_str() << '\n';
		if (boost::iequals(value, target))
			append(u);
	}
    this->set_up_geom_info();
}

void
GadgetronImagesVector::get_data(complex_float_t* data) const
{
	//std::copy(begin(), end(), data);
	//std::cout << "trying new const image wrap iterator...\n";
	GadgetronImagesVector::Iterator_const stop = end();
	GadgetronImagesVector::Iterator_const iter = begin();
	for (; iter != stop; ++iter, ++data)
		*data = (*iter).complex_float();
}

void
GadgetronImagesVector::set_data(const complex_float_t* data)
{
	//int dim[4];
	//size_t n = number();
	//get_image_dimensions(0, dim);
	//n *= dim[0];
	//n *= dim[1];
	//n *= dim[2];
	//n *= dim[3];
	//std::copy(data, data + n, begin());
	//std::cout << "trying new image wrap iterator...\n";
	GadgetronImagesVector::Iterator stop = end();
	GadgetronImagesVector::Iterator iter = begin();
	for (; iter != stop; ++iter, ++data)
		*iter = *data;
}

void
GadgetronImagesVector::get_real_data(float* data) const
{
	//std::cout << "in get_real_data...\n";
	GadgetronImagesVector::Iterator_const stop = end();
	GadgetronImagesVector::Iterator_const iter = begin();
	for (; iter != stop; ++iter, ++data)
		*data = *iter;
}

void
GadgetronImagesVector::set_real_data(const float* data)
{
	GadgetronImagesVector::Iterator stop = end();
	GadgetronImagesVector::Iterator iter = begin();
	for (; iter != stop; ++iter, ++data)
		*iter = *data;
}

static bool is_unit_vector(const float * const vec)
{
    return std::abs(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] - 1.F) < 1.e-4F;
}

static bool are_vectors_equal(const float * const vec1, const float * const vec2)
{
    for (int i=0; i<3; ++i)
        if (std::abs(vec1[i] - vec2[i]) > 1.e-4F)
            return false;
    return true;
}

static void print_slice_directions(const std::vector<gadgetron::shared_ptr<ImageWrap> > &images)
{
    std::cout << "\nGadgetronImagesVector::set_up_geom_info(): Slice direction alters between different slices. Expected it to be constant.\n";
    for (unsigned im=0; im<images.size(); ++im) {
        ISMRMRD::ImageHeader &ih = images[im]->head();
        float *sd = ih.slice_dir;
        std::cout << "Slice dir " << im << ": [" << sd[0] << ", " << sd[1] << ", " << sd[2] << "]\n";
    }
}

static void print_slice_distances(const std::vector<gadgetron::shared_ptr<ImageWrap> > &images)
{
    std::cout << "\nGadgetronImagesVector::set_up_geom_info(): Slice distances alters between slices. Expected it to be constant.\n";
    for (unsigned im=0; im<images.size()-1; ++im) {
        ISMRMRD::ImageHeader &ih1 = images[ im ]->head();
        ISMRMRD::ImageHeader &ih2 = images[im+1]->head();
        float projection_of_position_in_slice_dir_1 = ih1.position[0] * ih1.slice_dir[0] +
                ih1.position[1] * ih1.slice_dir[1] +
                ih1.position[2] * ih1.slice_dir[2];
        float projection_of_position_in_slice_dir_2 = ih2.position[0] * ih2.slice_dir[0] +
                ih2.position[1] * ih2.slice_dir[1] +
                ih2.position[2] * ih2.slice_dir[2];
        std::cout << "Spacing " << im << ": " << projection_of_position_in_slice_dir_1 - projection_of_position_in_slice_dir_2 << "\n";
    }
}

void
GadgetronImagesVector::print_header(const unsigned im_num)
{
    // Get image
    ISMRMRD::ImageHeader &ih = image_wrap(im_num).head();
    std::cout << "\n";
    std::cout << "phase:                  " << ih.phase                  << "\n";
    std::cout << "slice:                  " << ih.slice                  << "\n";
    std::cout << "average:                " << ih.average                << "\n";
    std::cout << "version:                " << ih.version                << "\n";
    std::cout << "channels:               " << ih.channels               << "\n";
    std::cout << "contrast:               " << ih.contrast               << "\n";
    std::cout << "data_type:              " << ih.data_type              << "\n";
    std::cout << "image_type:             " << ih.image_type             << "\n";
    std::cout << "repetition:             " << ih.repetition             << "\n";
    std::cout << "image_index:            " << ih.image_index            << "\n";
    std::cout << "measurement_uid:        " << ih.measurement_uid        << "\n";
    std::cout << "measurement_uid:        " << ih.measurement_uid        << "\n";
    std::cout << "image_series_index:     " << ih.image_series_index     << "\n";
    std::cout << "attribute_string_len:   " << ih.attribute_string_len   << "\n";
    std::cout << "acquisition_time_stamp: " << ih.acquisition_time_stamp << "\n";
    std::cout << "user_int:               "; for (int i=0;i<8;++i) std::cout << ih.user_int[i]               << " "; std::cout << "\n";
    std::cout << "user_float:             "; for (int i=0;i<8;++i) std::cout << ih.user_float[i]             << " "; std::cout << "\n";
    std::cout << "position:               "; for (int i=0;i<3;++i) std::cout << ih.position[i]               << " "; std::cout << "\n";
    std::cout << "read_dir:               "; for (int i=0;i<3;++i) std::cout << ih.read_dir[i]               << " "; std::cout << "\n";
    std::cout << "phase_dir:              "; for (int i=0;i<3;++i) std::cout << ih.phase_dir[i]              << " "; std::cout << "\n";
    std::cout << "slice_dir:              "; for (int i=0;i<3;++i) std::cout << ih.slice_dir[i]              << " "; std::cout << "\n";
    std::cout << "matrix_size:            "; for (int i=0;i<3;++i) std::cout << ih.matrix_size[i]            << " "; std::cout << "\n";
    std::cout << "field_of_view:          "; for (int i=0;i<3;++i) std::cout << ih.field_of_view[i]          << " "; std::cout << "\n";
    std::cout << "physiology_time_stamp:  "; for (int i=0;i<3;++i) std::cout << ih.physiology_time_stamp[i]  << " "; std::cout << "\n";
    std::cout << "patient_table_position: "; for (int i=0;i<3;++i) std::cout << ih.patient_table_position[i] << " "; std::cout << "\n";

    if (!acqs_info_.emtpty()) {
        std::cout << "XML data:\n";
        std::cout << acqs_info_.c_str() << "\n";
    }
}

void
GadgetronImagesVector::set_up_geom_info()
{
#ifndef NDEBUG
    std::cout << "\nSetting up geometrical info for GadgetronImagesVector...\n";
#endif

    if (number() < 1)
        return;

    if (!this->sorted())
        this->sort();

    // Get image
    ISMRMRD::ImageHeader &ih1 = image_wrap(0).head();

    // Check that the read, phase and slice directions are unit vectors and constant
    for (unsigned im=0; im<number(); ++im) {
        ISMRMRD::ImageHeader &ih = image_wrap(im).head();
        if (!(is_unit_vector(ih.read_dir) && is_unit_vector(ih.phase_dir) && is_unit_vector(ih.slice_dir))) {
            std::cout << "\nGadgetronImagesVector::set_up_geom_info(): read_dir, phase_dir and slice_dir should all be unit vectors.\n";
            return;
        }
        if (!(are_vectors_equal(ih1.read_dir,ih.read_dir) && are_vectors_equal(ih1.phase_dir,ih.phase_dir) && are_vectors_equal(ih1.slice_dir,ih.slice_dir))) {
            std::cout << "\nGadgetronImagesVector::set_up_geom_info(): read_dir, phase_dir and slice_dir should be constant over slices.\n";
            return;
        }
    }

    // Size
    // For the z-direction.
    // If it's a 3d image, matrix_size[2] == num voxels
    // If it's a 2d image, matrix_size[2] == 1, and number of slices is given by this->number()
    VoxelisedGeometricalInfo3D::Size size;
    for(unsigned i=0; i<3; ++i)
        size[i] = ih1.matrix_size[i];
    // If it's a stack of 2d images.
    if (size[2] == 1)
        size[2] = this->number();

    // The following will only work if the 0th index is read direction,
    // 1st is phase direction and 2nd is slice direction. This should be the case if sort has been called.

    // Spacing
    VoxelisedGeometricalInfo3D::Spacing spacing;
    for(unsigned i=0; i<3; ++i)
        spacing[i] = ih1.field_of_view[i] / size[i];

    // If there are more than 1 slices, then take the size of the voxel
    // in the z-direction to be the distance between voxel centres (this
    // accounts for under-sampled data (and also over-sampled).
    if (this->number() > 1) {

        // Calculate the spacing!
        ISMRMRD::ImageHeader &ih2 = image_wrap(1).head();
        float projection_of_position_in_slice_dir_1 = ih1.position[0] * ih1.slice_dir[0] +
                ih1.position[1] * ih1.slice_dir[1] +
                ih1.position[2] * ih1.slice_dir[2];
        float projection_of_position_in_slice_dir_2 = ih2.position[0] * ih2.slice_dir[0] +
                ih2.position[1] * ih2.slice_dir[1] +
                ih2.position[2] * ih2.slice_dir[2];
        spacing[2] = std::abs(projection_of_position_in_slice_dir_1 - projection_of_position_in_slice_dir_2);

        // Check: Loop over all images, and check that spacing is more-or-less constant
        for (unsigned im=0; im<number()-1; ++im) {

            ISMRMRD::ImageHeader &ih1 = image_wrap( im ).head();
            ISMRMRD::ImageHeader &ih2 = image_wrap(im+1).head();

            // 2. Check that spacing is constant
            float projection_of_position_in_slice_dir_1 = ih1.position[0] * ih1.slice_dir[0] +
                    ih1.position[1] * ih1.slice_dir[1] +
                    ih1.position[2] * ih1.slice_dir[2];
            float projection_of_position_in_slice_dir_2 = ih2.position[0] * ih2.slice_dir[0] +
                    ih2.position[1] * ih2.slice_dir[1] +
                    ih2.position[2] * ih2.slice_dir[2];
            float new_spacing = std::abs(projection_of_position_in_slice_dir_1 - projection_of_position_in_slice_dir_2);
            if (std::abs(spacing[2]-new_spacing) > 1.e-4F) {
                print_slice_distances(images_);
                return;
            }
        }
    }

    // Offset
    VoxelisedGeometricalInfo3D::Offset offset;
    for (int i=0; i<3; ++i)
        offset[i] = ih1.position[i];

    // Direction
    VoxelisedGeometricalInfo3D::DirectionMatrix direction;
    for (int axis=0; axis<3; ++axis) {
        direction[0][axis] = ih1.read_dir[axis];
        direction[1][axis] = ih1.phase_dir[axis];
        direction[2][axis] = ih1.slice_dir[axis];
    }

    // Initialise the geom info shared pointer
    _geom_info_sptr = std::make_shared<VoxelisedGeometricalInfo3D>
                (offset,spacing,size,direction);
}

void
CoilDataAsCFImage::get_data(float* re, float* im) const
{
	size_t n = img_.getNumberOfDataElements();
	const complex_float_t* ptr = img_.getDataPtr();
	for (size_t i = 0; i < n; i++) {
		complex_float_t z = ptr[i];
		re[i] = std::real(z);
		im[i] = std::imag(z);
	}
}

void 
CoilDataAsCFImage::set_data(const float* re, const float* im)
{
	size_t n = img_.getNumberOfDataElements();
	complex_float_t* ptr = img_.getDataPtr();
	for (size_t i = 0; i < n; i++)
		ptr[i] = complex_float_t((float)re[i], (float)im[i]);
}

void 
CoilImagesContainer::compute(MRAcquisitionData& ac)
{
	std::string par;
	ISMRMRD::IsmrmrdHeader header;
	ISMRMRD::Acquisition acq;
	par = ac.acquisitions_info();
	ISMRMRD::deserialize(par.c_str(), header);
	//ac.get_acquisition(0, acq);
	for (unsigned int i = 0; i < ac.number(); i++) {
		ac.get_acquisition(i, acq);
		if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_SLICE))
			break;
	}
	encoding_ = header.encoding[0];

	ISMRMRD::Encoding e = header.encoding[0];
	bool parallel = e.parallelImaging.is_present() &&
		e.parallelImaging().accelerationFactor.kspace_encoding_step_1 > 1;
	unsigned int nx = e.reconSpace.matrixSize.x;
	//unsigned int ny = e.reconSpace.matrixSize.y;
	//unsigned int nz = e.reconSpace.matrixSize.z;
	unsigned int ny = e.encodedSpace.matrixSize.y;
	unsigned int nz = e.encodedSpace.matrixSize.z;
	unsigned int nc = acq.active_channels();
	unsigned int readout = acq.number_of_samples();
	//std::cout << readout << '\n';
	//std::cout << nx << ' ' << ny << ' ' << nz << ' ' << nc << '\n';

	int nmap = 0;
	std::cout << "map ";

	for (unsigned int na = 0; na < ac.number();) {

		std::cout << ++nmap << ' ' << std::flush;

		std::vector<size_t> ci_dims;
		ci_dims.push_back(readout);
		ci_dims.push_back(ny);
		ci_dims.push_back(nz);
		ci_dims.push_back(nc);
		ISMRMRD::NDArray<complex_float_t> ci(ci_dims);
		memset(ci.getDataPtr(), 0, ci.getDataSize());

		int y = 0;
		for (;;) {
			ac.get_acquisition(na + y, acq);
			if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_SLICE))
				break;
			y++;
		}
		for (;;) {
			ac.get_acquisition(na + y, acq);
			int yy = acq.idx().kspace_encode_step_1;
			int zz = acq.idx().kspace_encode_step_2;
			//if (!e.parallelImaging.is_present() ||
			if (!parallel ||
				acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION) ||
				acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING)) {
				for (unsigned int c = 0; c < nc; c++) {
					for (unsigned int s = 0; s < readout; s++) {
						ci(s, yy, zz, c) = acq.data(s, c);
					}
				}
			}
			y++;
			if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE))
				break;
		}
		na += y;

		ifft3c(ci);

		shared_ptr<CoilData>
			sptr_ci(new CoilDataAsCFImage(readout, ny, nz, nc));
		CFImage& coil_im = (*(CoilDataAsCFImage*)sptr_ci.get()).image();
		memcpy(coil_im.getDataPtr(), ci.getDataPtr(), ci.getDataSize());
		append(sptr_ci);
	}
	std::cout << '\n';
}

void 
CoilSensitivitiesContainer::compute(CoilImagesContainer& cis)
{

	ISMRMRD::Encoding e = cis.encoding();
	unsigned int nx = e.reconSpace.matrixSize.x;
	//unsigned int ny = e.reconSpace.matrixSize.y;
	//unsigned int nz = e.reconSpace.matrixSize.z;
	unsigned int ny = e.encodedSpace.matrixSize.y;
	unsigned int nz = e.encodedSpace.matrixSize.z;
	int dim[4];
	cis(0).get_dim(dim);
	unsigned int readout = dim[0];
	unsigned int nc = dim[3];

	std::vector<size_t> cm_dims;
	cm_dims.push_back(readout);
	cm_dims.push_back(ny);
	cm_dims.push_back(nz);
	cm_dims.push_back(nc);
	ISMRMRD::NDArray<complex_float_t> cm(cm_dims);

	std::vector<size_t> csm_dims;
	csm_dims.push_back(nx);
	csm_dims.push_back(ny);
	csm_dims.push_back(nz);
	csm_dims.push_back(nc);
	ISMRMRD::NDArray<complex_float_t> csm(csm_dims);

	std::vector<size_t> img_dims;
	img_dims.push_back(nx);
	img_dims.push_back(ny);
	img_dims.push_back(nz);
	ISMRMRD::NDArray<float> img(img_dims);

	unsigned int nmap = 0;

	std::cout << "map ";
	for (nmap = 1; nmap <= cis.items(); nmap++) {
		std::cout << nmap << ' ' << std::flush;
		cis(nmap - 1).get_data(cm.getDataPtr());
		CoilData* ptr_img = new CoilDataAsCFImage(nx, ny, nz, nc);
		shared_ptr<CoilData> sptr_img(ptr_img);
		compute_csm_(cm, img, csm);
		ptr_img->set_data(csm.getDataPtr());
		append(sptr_img);
	}
	std::cout << '\n';
}

float 
CoilSensitivitiesContainer::max_(int nx, int ny, float* u)
{
	float r = 0.0;
	int i = 0;
	for (int iy = 0; iy < ny; iy++)
		for (int ix = 0; ix < nx; ix++, i++) {
			float t = fabs(u[i]);
			if (t > r)
				r = t;
		}
	return r;
}

void 
CoilSensitivitiesContainer::mask_noise_
(int nx, int ny, float* u, float noise, int* mask)
{
	int i = 0;
	for (int iy = 0; iy < ny; iy++)
		for (int ix = 0; ix < nx; ix++, i++) {
			float t = fabs(u[i]);
			mask[i] = (t > noise);
		}
}

int 
CoilSensitivitiesContainer::cleanup_mask_(int nx, int ny, int* mask, int bg, int minsz, int ex)
{
	int ll, il;
	int* listx = new int[nx*ny];
	int* listy = new int[nx*ny];
	int* inlist = new int[nx*ny];
	std::memset(inlist, 0, nx*ny * sizeof(int));
	int nc = 0;
	for (int iy = 0, i = 0; iy < ny; iy++) {
		for (int ix = 0; ix < nx; ix++, i++) {
			if (mask[i] == bg)
				continue;
			bool skip = false;
			ll = 1;
			listx[0] = ix;
			listy[0] = iy;
			inlist[i] = 1;
			il = 0;
			while (il < ll && ll < minsz) {
				int lx = listx[il];
				int ly = listy[il];
				int l = ll + ex;
				for (int jy = -l; jy <= l; jy++) {
					for (int jx = -l; jx <= l; jx++) {
						int kx = lx + jx;
						int ky = ly + jy;
						if (kx < 0 || kx >= nx)
							continue;
						if (ky < 0 || ky >= ny)
							continue;
						int j = kx + ky*nx;
						if (inlist[j])
							continue;
						if (mask[j] != bg) {
							listx[ll] = kx;
							listy[ll] = ky;
							inlist[j] = 1;
							ll++;
						}
					}
				}
				il++;
			}
			if (il == ll) {
				mask[i] = bg;
				nc++;
			}
			for (il = 0; il < ll; il++) {
				int lx = listx[il];
				int ly = listy[il];
				int j = lx + ly*nx;
				inlist[j] = 0;
			}
		}
	}
	//std::cout << nc << " mask pixels cleaned\n";
	delete[] listx;
	delete[] listy;
	delete[] inlist;
	return nc;
}

void 
CoilSensitivitiesContainer::smoothen_
(int nx, int ny, int nc,
	complex_float_t* u, complex_float_t* v,
	int* obj_mask, int w)
{
	const complex_float_t ONE(1.0, 0.0);
	const complex_float_t TWO(2.0, 0.0);
	for (int ic = 0, i = 0; ic < nc; ic++)
		for (int iy = 0, k = 0; iy < ny; iy++)
			for (int ix = 0; ix < nx; ix++, i++, k++) {
				//if (edge_mask[k]) {
				//	v[i] = u[i];
				//	continue;
				//}
				//if (!obj_mask[k]) {
				//	v[i] = 0;
				//	continue;
				//}
				int n = 0;
				complex_float_t r(0.0, 0.0);
				complex_float_t s(0.0, 0.0);
				for (int jy = -w; jy <= w; jy++)
					for (int jx = -w; jx <= w; jx++) {
						if (ix + jx < 0 || ix + jx >= nx)
							continue;
						if (iy + jy < 0 || iy + jy >= ny)
							continue;
						int j = i + jx + jy*nx;
						int l = k + jx + jy*nx;
						if (i != j && obj_mask[l]) { // == obj_mask[k]) { // && !edge_mask[l]) {
							n++;
							r += ONE;
							s += u[j];
						}
					}
				//v[i] = obj_mask[k];
				if (n > 0)
					v[i] = (u[i] + s / r) / TWO;
				else
					v[i] = u[i];
			}
	memcpy(u, v, nx*ny*nc * sizeof(complex_float_t));
}

void 
CoilSensitivitiesContainer::compute_csm_(
	ISMRMRD::NDArray<complex_float_t>& cm,
	ISMRMRD::NDArray<float>& img,
	ISMRMRD::NDArray<complex_float_t>& csm
)
{
	int ndims = cm.getNDim();
	const size_t* dims = cm.getDims();
	unsigned int readout = (unsigned int)dims[0];
	unsigned int ny = (unsigned int)dims[1];
	unsigned int nz = (unsigned int)dims[2];
	unsigned int nc = (unsigned int)dims[3];
	unsigned int nx = (unsigned int)img.getDims()[0];

	std::vector<size_t> cm0_dims;
	cm0_dims.push_back(nx);
	cm0_dims.push_back(ny);
	cm0_dims.push_back(nz);
	cm0_dims.push_back(nc);

	ISMRMRD::NDArray<complex_float_t> cm0(cm0_dims);
	for (unsigned int c = 0; c < nc; c++) {
		for (unsigned int z = 0; z < nz; z++) {
			for (unsigned int y = 0; y < ny; y++) {
				for (unsigned int x = 0; x < nx; x++) {
					uint16_t xout = x + (readout - nx) / 2;
					cm0(x, y, z, c) = cm(xout, y, z, c);
				}
			}
		}
	}

	int* object_mask = new int[nx*ny*nc];
	memset(object_mask, 0, nx*ny*nc * sizeof(int));

	ISMRMRD::NDArray<complex_float_t> w(cm0);

	float* ptr_img = img.getDataPtr();
	for (unsigned int z = 0; z < nz; z++) {
		for (unsigned int y = 0; y < ny; y++) {
			for (unsigned int x = 0; x < nx; x++) {
				float r = 0.0;
				for (unsigned int c = 0; c < nc; c++) {
					float s = std::abs(cm0(x, y, z, c));
					r += s*s;
				}
				img(x, y, z) = (float)std::sqrt(r);
			}
		}
	}

	float max_im = max_(nx, ny, ptr_img);
	float noise = max_(5, 5, ptr_img) + (float)1e-6*max_im;
	//std::cout << "max_im: " << max_im << ", noise: " << noise << '\n';
	mask_noise_(nx, ny, ptr_img, noise, object_mask);
	//for (int i = 2; i < 20; i++) {
	//	int cleaned = cleanup_mask_(nx, ny, object_mask, 0, i, 0);
	//	if (cleaned < 1)
	//		break;
	//}
	//for (int i = 2; i < 20; i++) {
	//	int cleaned = cleanup_mask_(nx, ny, object_mask, 1, i, 0);
	//	if (cleaned < 1)
	//		break;
	//}
	cleanup_mask_(nx, ny, object_mask, 0, 2, 0);
	cleanup_mask_(nx, ny, object_mask, 0, 3, 0);
	cleanup_mask_(nx, ny, object_mask, 0, 4, 0);
	//cleanup_mask_(nx, ny, object_mask, 1, 2, 0);
	//cleanup_mask_(nx, ny, object_mask, 1, 3, 0);
	//cleanup_mask_(nx, ny, object_mask, 1, 4, 0);

	for (int i = 0; i < csm_smoothness_; i++)
		smoothen_(nx, ny, nc, cm0.getDataPtr(), w.getDataPtr(), object_mask, 1);

	for (unsigned int z = 0; z < nz; z++) {
		for (unsigned int y = 0; y < ny; y++) {
			for (unsigned int x = 0; x < nx; x++) {
				float r = 0.0;
				for (unsigned int c = 0; c < nc; c++) {
					float s = std::abs(cm0(x, y, z, c));
					r += s*s;
				}
				img(x, y, z) = (float)std::sqrt(r);
			}
		}
	}

	for (unsigned int z = 0; z < nz; z++) {
		for (unsigned int y = 0; y < ny; y++) {
			for (unsigned int x = 0; x < nx; x++) {
				float r = img(x, y, z);
				float s;
				if (r != 0.0)
					s = (float)(1.0 / r);
				else
					s = 0.0;
				complex_float_t zs(s, 0.0);
				for (unsigned int c = 0; c < nc; c++) {
					csm(x, y, z, c) = cm0(x, y, z, c) * zs;
				}
			}
		}
	}

	delete[] object_mask;

}

CoilSensitivitiesAsImages::CoilSensitivitiesAsImages(const char* file)
{
	Mutex mtx;
	mtx.lock();
	ISMRMRD::Dataset csm_file(file, "dataset");
	int nm = csm_file.getNumberOfImages("csm");
	mtx.unlock();
	for (int i = 0; i < nm; i++) {
		shared_ptr<CoilData> sptr_img(new CoilDataAsCFImage);
		mtx.lock();
		CFImage& csm = (*(CoilDataAsCFImage*)sptr_img.get()).image();
		csm_file.readImage("csm", i, csm);
		mtx.unlock();
		append(sptr_img);
	}
	csm_smoothness_ = 0;
}
