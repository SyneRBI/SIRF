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
\ingroup Gadgetron Data Containers
\brief Implementation file for SIRF data container classes for Gadgetron data.

\author Evgueni Ovtchinnikov
\author CCP PETMR
*/

#include "gadgetron_data_containers.h"

void 
ImageWrap::get_cmplx_data(double* re, double* im) const
{
	int dim[4];
	size_t n = get_dim(dim);
	if (type_ == ISMRMRD::ISMRMRD_CXFLOAT) {
		const CFImage& img = *(const CFImage*)ptr_;
		const complex_float_t* ptr = img.getDataPtr();
		for (size_t i = 0; i < n; i++) {
			complex_float_t z = ptr[i];
			re[i] = std::real(z);
			im[i] = std::imag(z);
		}
	}
	else if (type_ == ISMRMRD::ISMRMRD_CXDOUBLE) {
		const CDImage& img = *(const CDImage*)ptr_;
		const complex_double_t* ptr = img.getDataPtr();
		for (size_t i = 0; i < n; i++) {
			complex_double_t z = ptr[i];
			re[i] = std::real(z);
			im[i] = std::imag(z);
		}
	}
	else {
		get_data(re);
		for (size_t i = 0; i < n; i++)
			im[i] = 0;
	}
}

void 
ImageWrap::set_cmplx_data(const double* re, const double* im) const
{
	int dim[4];
	size_t n = get_dim(dim);
	if (type_ == ISMRMRD::ISMRMRD_CXFLOAT) {
		CFImage& img = *(CFImage*)ptr_;
		complex_float_t* ptr = img.getDataPtr();
		for (size_t i = 0; i < n; i++)
			ptr[i] = std::complex<float>((float)re[i], (float)im[i]);
	}
	else if (type_ == ISMRMRD::ISMRMRD_CXDOUBLE) {
		CDImage& img = *(CDImage*)ptr_;
		complex_double_t* ptr = img.getDataPtr();
		for (size_t i = 0; i < n; i++)
			ptr[i] = std::complex<double>(re[i], im[i]);
	}
}

int 
AcquisitionsContainer::get_acquisitions_dimensions(size_t ptr_dim)
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
	int nrd = ordered() ? 3 : 2;
	// number of regular readouts
	int nrr = 0;
	//int not_reg = 0;
	for (; y < na;) {
		for (; y < na && ordered();) {
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
			if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE) && ordered())
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
	//dim[0] = acq.number_of_samples();
	//dim[1] = acq.active_channels();
	//dim[2] = my; // e.reconSpace.matrixSize.y;
	//dim[3] = slice;
	//return not_reg;
}

int 
AcquisitionsContainer::set_acquisitions_data
(boost::shared_ptr<AcquisitionsContainer> sptr_ac,
	int na, int nc, int ns, const double* re, const double* im)
{
	sptr_ac->set_parameters(par_);
	sptr_ac->write_parameters();
	sptr_ac->ordered_ = ordered();
	ISMRMRD::Acquisition acq;
	int ma = number();
	for (int a = 0, i = 0; a < ma; a++) {
		get_acquisition(a, acq);
		if (TO_BE_IGNORED(acq) && ma > na) {
			std::cout << "ignoring acquisition " << a << '\n';
			continue;
		}
		unsigned int mc = acq.active_channels();
		unsigned int ms = acq.number_of_samples();
		if (mc != nc || ms != ns)
			return -1;
		//for (size_t c = 0; c < nc; c++) {
		//	for (size_t s = 0; s < ns; s++, i++) {
		//		acq.data(s, c) = complex_float_t((float)re[i], (float)im[i]);
		//	}
		//}
		for (int c = 0; c < nc; c++)
			for (int s = 0; s < ns; s++, i++)
				acq.data(s, c) = complex_float_t((float)re[i], (float)im[i]);
		sptr_ac->append_acquisition(acq);
	}
	return 0;
}

void 
AcquisitionsContainer::get_acquisitions_flags(unsigned int n, int* flags)
{
	ISMRMRD::Acquisition acq;
	unsigned int na = number();
	for (unsigned int a = 0, i = 0; a < na; a++) {
		get_acquisition(a, acq);
		if (TO_BE_IGNORED(acq) && n < na) {
			std::cout << "ignoring acquisition " << a << '\n';
			continue;
		}
		flags[i++] = (int)acq.flags();
	}
}

unsigned int 
AcquisitionsContainer::get_acquisitions_data(unsigned int slice, double* re, double* im)
{
	ISMRMRD::Acquisition acq;
	unsigned int na = number();
	unsigned int n = 0;
	if (slice >= na) {
		for (unsigned int a = 0, i = 0; a < na; a++) {
			get_acquisition(a, acq);
			if (TO_BE_IGNORED(acq) && slice > na) {
				std::cout << "ignoring acquisition " << a << '\n';
				continue;
			}
			n++;
			unsigned int nc = acq.active_channels();
			unsigned int ns = acq.number_of_samples();
			for (unsigned int c = 0; c < nc; c++) {
				for (unsigned int s = 0; s < ns; s++, i++) {
					complex_float_t z = acq.data(s, c);
					re[i] = std::real(z);
					im[i] = std::imag(z);
				}
			}
		}
		return n;
	}
	int* dim = new int[3];
	size_t ptr_dim = (size_t)dim;
	get_acquisitions_dimensions(ptr_dim);
	unsigned int ny = dim[2]; //e.reconSpace.matrixSize.y;
							  //unsigned int ny = dim[1]; //e.reconSpace.matrixSize.y;
	delete[] dim;
	unsigned int y = 0;
	for (; y + ny*slice < na;) {
		get_acquisition(y + ny*slice, acq);
		if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_SLICE))
			break;
		y++;
	}
	for (; y + ny*slice < na; n++) {
		get_acquisition(y + ny*slice, acq);
		unsigned int nc = acq.active_channels();
		unsigned int ns = acq.number_of_samples();
		for (unsigned int c = 0; c < nc; c++) {
			for (unsigned int s = 0; s < ns; s++) {
				complex_float_t z = acq.data(s, c);
				re[s + ns*(n + ny*c)] = std::real(z);
				im[s + ns*(n + ny*c)] = std::imag(z);
			}
		}
		y++;
		if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE))
			break;
	}
	return n;
}

void 
AcquisitionsContainer::axpby
(complex_double_t a, const ISMRMRD::Acquisition& acq_x,
	complex_double_t b, ISMRMRD::Acquisition& acq_y)
{
	complex_float_t* px;
	complex_float_t* py;
	for (px = acq_x.data_begin(), py = acq_y.data_begin();
		px != acq_x.data_end() && py != acq_y.data_end(); px++, py++) {
		if (b == complex_double_t(0.0))
			*py = a*complex_double_t(*px);
		else
			*py = a*complex_double_t(*px) + b*complex_double_t(*py);
	}
}
