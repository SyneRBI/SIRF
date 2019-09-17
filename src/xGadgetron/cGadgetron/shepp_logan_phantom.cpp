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

Contains code snippets from ISMRMRD/utilities/generate_cartesian_shepp_logan.cpp
by Michael S. Hansen

ISMRMRD SOFTWARE LICENSE JULY 2013

PERMISSION IS HEREBY GRANTED, FREE OF CHARGE, TO ANY PERSON OBTAINING
A COPY OF THIS SOFTWARE AND ASSOCIATED DOCUMENTATION FILES (THE
"SOFTWARE"), TO DEAL IN THE SOFTWARE WITHOUT RESTRICTION, INCLUDING
WITHOUT LIMITATION THE RIGHTS TO USE, COPY, MODIFY, MERGE, PUBLISH,
DISTRIBUTE, SUBLICENSE, AND/OR SELL COPIES OF THE SOFTWARE, AND TO
PERMIT PERSONS TO WHOM THE SOFTWARE IS FURNISHED TO DO SO, SUBJECT TO
THE FOLLOWING CONDITIONS:

THE ABOVE COPYRIGHT NOTICE, THIS PERMISSION NOTICE, AND THE LIMITATION
OF LIABILITY BELOW SHALL BE INCLUDED IN ALL COPIES OR REDISTRIBUTIONS
OF SUBSTANTIAL PORTIONS OF THE SOFTWARE.

SOFTWARE IS BEING DEVELOPED IN PART AT THE NATIONAL HEART, LUNG, AND BLOOD
INSTITUTE, NATIONAL INSTITUTES OF HEALTH BY AN EMPLOYEE OF THE FEDERAL
GOVERNMENT IN THE COURSE OF HIS OFFICIAL DUTIES. PURSUANT TO TITLE 17,
SECTION 105 OF THE UNITED STATES CODE, THIS SOFTWARE IS NOT SUBJECT TO
COPYRIGHT PROTECTION AND IS IN THE PUBLIC DOMAIN. EXCEPT AS CONTAINED IN
THIS NOTICE, THE NAME OF THE AUTHORS, THE NATIONAL HEART, LUNG, AND BLOOD
INSTITUTE (NHLBI), OR THE NATIONAL INSTITUTES OF HEALTH (NIH) MAY NOT
BE USED TO ENDORSE OR PROMOTE PRODUCTS DERIVED FROM THIS SOFTWARE WITHOUT
SPECIFIC PRIOR WRITTEN PERMISSION FROM THE NHLBI OR THE NIH.THE SOFTWARE IS
PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


*/

#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/xml.h"
#include "ismrmrd/dataset.h"
#include "ismrmrd/version.h"
#include "sirf/Gadgetron/ismrmrd_phantom.h"
#include "sirf/Gadgetron/ismrmrd_fftw.h"
#include "sirf/Gadgetron/cgadgetron_shared_ptr.h"

using namespace ISMRMRD;
using namespace gadgetron;

void generate_cartesian_shepp_logan(unsigned int matrix_size, unsigned int ncoils, unsigned int ros, std::string file)
{
	shared_ptr<NDArray<complex_float_t> > phantom = shepp_logan_phantom(matrix_size);
	shared_ptr<NDArray<complex_float_t> > coils = generate_birdcage_sensititivies(matrix_size, ncoils, 1.5);

	std::vector<size_t> dims;
	dims.push_back(matrix_size*ros); //oversampling in the readout direction
	dims.push_back(matrix_size);
	dims.push_back(ncoils);

	NDArray<complex_float_t> coil_images(dims);
	memset(coil_images.getDataPtr(), 0, coil_images.getDataSize());

	size_t readout = matrix_size*ros;

	for (unsigned int c = 0; c < ncoils; c++) {
		for (unsigned int y = 0; y < matrix_size; y++) {
			for (unsigned int x = 0; x < matrix_size; x++) {
				uint16_t xout = x + (readout - matrix_size) / 2;
				coil_images(xout, y, c) = (*phantom)(x, y) * (*coils)(x, y, c);
			}
		}
	}

	Dataset d(file.c_str(), "dataset", true);
	Acquisition acq;
	memset((void*)acq.getDataPtr(), 0, acq.getDataSize());
	acq.available_channels() = ncoils;
	acq.center_sample() = (readout >> 1);

}