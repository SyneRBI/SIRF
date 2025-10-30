#========================================================================
# Copyright 2021 Science Technology Facilities Council
#
# This file is part of the SyneRBI Synergistic Image Reconstruction Framework (SIRF).
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#         http://www.apache.org/licenses/LICENSE-2.0.txt
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#=========================================================================

import os, shutil, numpy
import unittest
import sirf.STIR as pet
from sirf.Utilities import examples_data_path

pet.AcquisitionData.set_storage_scheme('file')
pet.set_verbosity(0)


class TestSTIRObjectiveFunction(unittest.TestCase):

    def setUp(self):
        data_path = os.path.join(examples_data_path('PET'), 'thorax_single_slice')

        image = pet.ImageData(os.path.join(data_path,'emission.hv'))

        am = pet.AcquisitionModelUsingRayTracingMatrix()
        am.set_num_tangential_LORs(5)
        templ = pet.AcquisitionData(os.path.join(data_path,'template_sinogram.hs'))
        am.set_up(templ,image)
        acquired_data=am.forward(image)
        am.set_background_term(acquired_data*0 + numpy.mean(acquired_data.as_array()))
        am.set_up(templ,image)
        obj_fun = pet.make_Poisson_loglikelihood(acquired_data)
        obj_fun.set_acquisition_model(am)
        prior = pet.QuadraticPrior()
        prior.set_penalisation_factor(20)
        obj_fun.set_prior(prior)
        obj_fun.set_up(image)

        self.obj_fun = obj_fun
        self.image = image

    #def tearDown(self):


    def test_Poisson_loglikelihood_call(self):
        x = self.image * 0 + 1
        a = self.obj_fun.value(x)
        b = self.obj_fun(x)

        numpy.testing.assert_almost_equal(a,b)

    def test_Hessian(self, subset=-1, eps=1e-3):
        """Checks that grad(x + dx) - grad(x) is close to H(x)*dx
        """
        x = self.image
        dx = x.clone()
        dx *= eps
        dx += eps/2
        y = x + dx
        gx = self.obj_fun.gradient(x, subset)
        gy = self.obj_fun.gradient(y, subset)
        dg = gy - gx
        Hdx = self.obj_fun.multiply_with_Hessian(x, dx, subset)
        q = (dg - Hdx).norm()/dg.norm()
        print('norm of (x): %f' % x.norm())
        print('norm of (x + dx): %f' % y.norm())
        print('norm of grad(x): %f' % gx.norm())
        print('norm of grad(x + dx): %f' % gy.norm())
        print('norm of grad(x + dx) - grad(x): %f' % dg.norm())
        print('norm of H(x)*dx: %f' % Hdx.norm())
        print('relative difference: %f' % q)
        assert q <= .002

    def test_gradient_inplace(self):
        '''Checks that gradient computed in-place returns the same as without in-place.

        Note: the absolute tolerance of the test was increased from the default 1e-8 to 1e-5, see 'test_gradient_out'
        '''
        x = self.image
        g1 = self.obj_fun.gradient(x)
        self.obj_fun.gradient(x, out=x)
        numpy.testing.assert_allclose(x.asarray(), g1.asarray(), atol=1e-5)

    def test_gradient_out(self):
        '''Checks that gradient with 'out' parameter returns the same as without the 'out' parameter.

        Note: the absolute tolerance of the test was increased from the default 1e-8 to 1e-5.
        https://github.com/SyneRBI/SIRF/issues/1349
        '''
        x = self.image
        g1 = self.obj_fun.gradient(x)
        g2 = g1 * 0
        self.obj_fun.gradient(x, out=g2)
        numpy.testing.assert_allclose(g2.asarray(), g1.asarray(), atol=1e-5)
