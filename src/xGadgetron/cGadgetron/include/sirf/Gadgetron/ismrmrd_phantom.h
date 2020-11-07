/*
 * ismrmrd_phantom.h

 *  Created on: Apr 1, 2013
 *      Author: Michael S. Hansen

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

#include <boost/shared_ptr.hpp>
#include <vector>
#include <complex>
#include "ismrmrd/ismrmrd.h"
#include "sirf/Gadgetron/cgadgetron_shared_ptr.h"

#ifndef ISMRMRD_PHANTOM_H_
#define ISMRMRD_PHANTOM_H_

namespace ISMRMRD {

	class PhantomEllipse
	{
	public:
		PhantomEllipse(float A, float a, float b, float x0, float y0, float phi)
			: A_(A)
			, a_(a)
			, b_(b)
			, x0_(x0)
			, y0_(y0)
			, phi_(phi)
		{

		}

		bool isInside(float x, float y)
		{
			float asq = a_*a_;                // a^2
			float bsq = b_*b_;                // b^2
			float phi = phi_*3.14159265359 / 180; // rotation angle in radians
			float x0 = x - x0_; 					  // x offset
			float y0 = y - y0_;                     // y offset
			float cosp = cos(phi);
			float sinp = sin(phi);
			return (((x0*cosp + y0*sinp)*(x0*cosp + y0*sinp)) / asq + ((y0*cosp - x0*sinp)*(y0*cosp - x0*sinp)) / bsq <= 1);
		}

		float getAmplitude()
		{
			return A_;
		}


	protected:

		float A_;
		float a_;
		float b_;
		float x0_;
		float y0_;
		float phi_;
	};

	gadgetron::shared_ptr< std::vector<PhantomEllipse> > shepp_logan_ellipses();
	gadgetron::shared_ptr< std::vector<PhantomEllipse> > modified_shepp_logan_ellipses();
	gadgetron::shared_ptr<NDArray<complex_float_t> > phantom(std::vector<PhantomEllipse>& coefficients, unsigned int matrix_size);
	gadgetron::shared_ptr<NDArray<complex_float_t> > shepp_logan_phantom(unsigned int matrix_size);
	gadgetron::shared_ptr<NDArray<complex_float_t> > generate_birdcage_sensititivies(unsigned int matrix_size, unsigned int ncoils, float relative_radius);
	int add_noise(NDArray<complex_float_t> & a, float sd);
	int add_noise(Acquisition & a, float sd);

};

#endif /* ISMRMRD_PHANTOM_H_ */
