#========================================================================
# Copyright 2017 - 2021 University College London
# Copyright 2017 - 2021 Science Technology Facilities Council
#
# This file is part of the SyneRBI Synergistic Image Reconstruction Framework (SIRF) SuperBuild.
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

import numpy
import os
import shutil
import unittest
import sirf.STIR as pet
import sirf.Gadgetron as mr
from sirf.Utilities import  examples_data_path


class TestSTIRImageDataAlgebra(unittest.TestCase):

    def setUp(self):
        os.chdir(examples_data_path('PET'))
        # Copy files to a working folder and change directory to where these files are.
        # We do this to avoid cluttering your SIRF files. This way, you can delete
        # working_folder and start from scratch.
        shutil.rmtree('working_folder/brain',True)
        shutil.copytree('brain','working_folder/brain')
        os.chdir('working_folder/brain')

        self.cwd = os.getcwd()

    
    def tearDown(self):
        shutil.rmtree(self.cwd)

    def test_divide_scalar(self):
        os.chdir(self.cwd)
        image1 = pet.ImageData('emission.hv')
        image2 = pet.ImageData('emission.hv')
        image1.fill(1.)
        image2.fill(2.)
        
        tmp = image1/1.
        numpy.testing.assert_array_equal(image1.as_array(), tmp.as_array())
    
        tmp1 = image1.divide(1.)
        numpy.testing.assert_array_equal(tmp.as_array(), tmp1.as_array())
        
        image1.divide(1., out=image2)
        numpy.testing.assert_array_equal(tmp.as_array(), image2.as_array())
        
    def test_divide_datacontainer(self):
        os.chdir(self.cwd)
        # add 1 because the data contains zeros and divide is not going to be happy
        image1 = pet.ImageData('emission.hv') + 1
        image2 = pet.ImageData('emission.hv') + 1
                
        tmp = image1/image2

        numpy.testing.assert_array_almost_equal(
            numpy.ones(image1.shape, dtype=numpy.float32), tmp.as_array()
            )
    
        tmp1 = image1.divide(image2)
        numpy.testing.assert_array_almost_equal(
            numpy.ones(image1.shape, dtype=numpy.float32), tmp1.as_array()
            )
        
        tmp1.fill(2.)
        image1.divide(image2, out=tmp1)
        
        numpy.testing.assert_array_almost_equal(
            numpy.ones(image1.shape, dtype=numpy.float32), tmp1.as_array()
            )
        

    def test_multiply_scalar(self):
        os.chdir(self.cwd)
        image1 = pet.ImageData('emission.hv')
        image2 = pet.ImageData('emission.hv')
        
        image2.fill(2.)
        
        tmp = image1 * 1.
        numpy.testing.assert_array_equal(image1.as_array(), tmp.as_array())
    
        tmp1 = image1.multiply(1.)
        numpy.testing.assert_array_equal(tmp.as_array(), tmp1.as_array())
        
        image1.multiply(1., out=image2)
        numpy.testing.assert_array_equal(tmp.as_array(), image2.as_array())
        
    def test_multiply_datacontainer(self):
        os.chdir(self.cwd)
        image1 = pet.ImageData('emission.hv')
        image2 = pet.ImageData('emission.hv')
        
        image2.fill(1.)
        tmp = image1 * image2

        numpy.testing.assert_array_almost_equal(
            image1.as_array(), tmp.as_array()
            )
    
        tmp1 = image1.multiply(image2)
        numpy.testing.assert_array_almost_equal(
            image1.as_array(), tmp1.as_array()
            )
        
        tmp1.fill(2.)
        image1.multiply(image2, out=tmp1)
        
        numpy.testing.assert_array_almost_equal(
            image1.as_array(), tmp1.as_array()
            )
        
    def test_add_scalar(self):
        os.chdir(self.cwd)
        image1 = pet.ImageData('emission.hv')
        image2 = pet.ImageData('emission.hv')
        image1.fill(0)
        image2.fill(1)
        
        tmp = image1 + 1.
        numpy.testing.assert_array_equal(image2.as_array(), tmp.as_array())
    
        tmp1 = image1.add(1.)
        numpy.testing.assert_array_equal(tmp.as_array(), tmp1.as_array())
        
        tmp1.fill(0)
        image1.add(1., out=tmp1)
        numpy.testing.assert_array_equal(tmp1.as_array(), image2.as_array())
    
    def test_add_datacontainer(self):
        os.chdir(self.cwd)
        image1 = pet.ImageData('emission.hv')
        image2 = pet.ImageData('emission.hv')
        image1.fill(0.)
        image2.fill(1.)
        tmp = image1 + image2

        numpy.testing.assert_array_almost_equal(
            numpy.ones(image1.shape, dtype=numpy.float32), tmp.as_array()
            )
    
        tmp1 = image1.add(image2)
        
        numpy.testing.assert_array_almost_equal(
            numpy.ones(image1.shape, dtype=numpy.float32), tmp1.as_array()
            )
        
        tmp1.fill(2.)
        image1.add(image2, out=tmp1)
        
        numpy.testing.assert_array_almost_equal(
            numpy.ones(image1.shape, dtype=numpy.float32), tmp1.as_array()
            )
        
    
    def test_subtract_scalar(self):
        os.chdir(self.cwd)
        image1 = pet.ImageData('emission.hv')
        image2 = pet.ImageData('emission.hv')
        image1.fill(2)
        image2.fill(1)
        
        tmp = image1 - 1.
        numpy.testing.assert_array_equal(image2.as_array(), tmp.as_array())
    
        tmp1 = image1.subtract(1.)
        numpy.testing.assert_array_equal(tmp.as_array(), tmp1.as_array())
        
        tmp1.fill(0)
        image1.subtract(1., out=tmp1)
        numpy.testing.assert_array_equal(tmp1.as_array(), image2.as_array())
    
    def test_subtract_datacontainer(self):
        os.chdir(self.cwd)
        image1 = pet.ImageData('emission.hv')
        image2 = pet.ImageData('emission.hv')
        
        tmp = image1 - image2

        numpy.testing.assert_array_almost_equal(
            numpy.zeros(image1.shape, dtype=numpy.float32), tmp.as_array()
            )
    
        tmp1 = image1.subtract(image2)
        
        numpy.testing.assert_array_almost_equal(
            numpy.zeros(image1.shape, dtype=numpy.float32), tmp1.as_array()
            )
        
        tmp1.fill(2.)
        image1.subtract(image2, out=tmp1)
        
        numpy.testing.assert_array_almost_equal(
            numpy.zeros(image1.shape, dtype=numpy.float32), tmp1.as_array()
            )

class TestSTIRAcquisitionDataAlgebra(unittest.TestCase):

    def setUp(self):
        os.chdir(examples_data_path('PET'))
        # Copy files to a working folder and change directory to where these files are.
        # We do this to avoid cluttering your SIRF files. This way, you can delete
        # working_folder and start from scratch.
        shutil.rmtree('working_folder/brain',True)
        shutil.copytree('brain','working_folder/brain')
        os.chdir('working_folder/brain')

        self.cwd = os.getcwd()

    
    def tearDown(self):
        shutil.rmtree(self.cwd)

    def test_divide_scalar(self):
        os.chdir(self.cwd)
        image1 = pet.ImageData('attenuation.hv')
        image2 = pet.ImageData('attenuation.hv')
        image1.fill(1.)
        image2.fill(2.)
        
        tmp = image1/1.
        numpy.testing.assert_array_equal(image1.as_array(), tmp.as_array())
    
        tmp1 = image1.divide(1.)
        numpy.testing.assert_array_equal(tmp.as_array(), tmp1.as_array())
        
        image1.divide(1., out=image2)
        numpy.testing.assert_array_equal(tmp.as_array(), image2.as_array())
        
    def test_divide_datacontainer(self):
        os.chdir(self.cwd)
        # add 1 because the data contains zeros and divide is not going to be happy
        image1 = pet.ImageData('attenuation.hv') + 1
        image2 = pet.ImageData('attenuation.hv') + 1
                
        tmp = image1/image2

        numpy.testing.assert_array_almost_equal(
            numpy.ones(image1.shape, dtype=numpy.float32), tmp.as_array()
            )
    
        tmp1 = image1.divide(image2)
        numpy.testing.assert_array_almost_equal(
            numpy.ones(image1.shape, dtype=numpy.float32), tmp1.as_array()
            )
        
        tmp1.fill(2.)
        image1.divide(image2, out=tmp1)
        
        numpy.testing.assert_array_almost_equal(
            numpy.ones(image1.shape, dtype=numpy.float32), tmp1.as_array()
            )
        

    def test_multiply_scalar(self):
        os.chdir(self.cwd)
        image1 = pet.ImageData('attenuation.hv')
        image2 = pet.ImageData('attenuation.hv')
        
        image2.fill(2.)
        
        tmp = image1 * 1.
        numpy.testing.assert_array_equal(image1.as_array(), tmp.as_array())
    
        tmp1 = image1.multiply(1.)
        numpy.testing.assert_array_equal(tmp.as_array(), tmp1.as_array())
        
        image1.multiply(1., out=image2)
        numpy.testing.assert_array_equal(tmp.as_array(), image2.as_array())
        
    def test_multiply_datacontainer(self):
        os.chdir(self.cwd)
        image1 = pet.ImageData('attenuation.hv')
        image2 = pet.ImageData('attenuation.hv')
        
        image2.fill(1.)
        tmp = image1 * image2

        numpy.testing.assert_array_almost_equal(
            image1.as_array(), tmp.as_array()
            )
    
        tmp1 = image1.multiply(image2)
        numpy.testing.assert_array_almost_equal(
            image1.as_array(), tmp1.as_array()
            )
        
        tmp1.fill(2.)
        image1.multiply(image2, out=tmp1)
        
        numpy.testing.assert_array_almost_equal(
            image1.as_array(), tmp1.as_array()
            )
        
    def test_add_scalar(self):
        os.chdir(self.cwd)
        image1 = pet.ImageData('attenuation.hv')
        image2 = pet.ImageData('attenuation.hv')
        image1.fill(0)
        image2.fill(1)
        
        tmp = image1 + 1.
        numpy.testing.assert_array_equal(image2.as_array(), tmp.as_array())
    
        tmp1 = image1.add(1.)
        numpy.testing.assert_array_equal(tmp.as_array(), tmp1.as_array())
        
        tmp1.fill(0)
        image1.add(1., out=tmp1)
        numpy.testing.assert_array_equal(tmp1.as_array(), image2.as_array())
    
    def test_add_datacontainer(self):
        os.chdir(self.cwd)
        image1 = pet.ImageData('attenuation.hv')
        image2 = pet.ImageData('attenuation.hv')
        image1.fill(0.)
        image2.fill(1.)
        tmp = image1 + image2

        numpy.testing.assert_array_almost_equal(
            numpy.ones(image1.shape, dtype=numpy.float32), tmp.as_array()
            )
    
        tmp1 = image1.add(image2)
        
        numpy.testing.assert_array_almost_equal(
            numpy.ones(image1.shape, dtype=numpy.float32), tmp1.as_array()
            )
        
        tmp1.fill(2.)
        image1.add(image2, out=tmp1)
        
        numpy.testing.assert_array_almost_equal(
            numpy.ones(image1.shape, dtype=numpy.float32), tmp1.as_array()
            )
        
    
    def test_subtract_scalar(self):
        os.chdir(self.cwd)
        image1 = pet.ImageData('attenuation.hv')
        image2 = pet.ImageData('attenuation.hv')
        image1.fill(2)
        image2.fill(1)
        
        tmp = image1 - 1.
        numpy.testing.assert_array_equal(image2.as_array(), tmp.as_array())
    
        tmp1 = image1.subtract(1.)
        numpy.testing.assert_array_equal(tmp.as_array(), tmp1.as_array())
        
        tmp1.fill(0)
        image1.subtract(1., out=tmp1)
        numpy.testing.assert_array_equal(tmp1.as_array(), image2.as_array())
    
    def test_subtract_datacontainer(self):
        os.chdir(self.cwd)
        image1 = pet.ImageData('attenuation.hv')
        image2 = pet.ImageData('attenuation.hv')
        
        tmp = image1 - image2

        numpy.testing.assert_array_almost_equal(
            numpy.zeros(image1.shape, dtype=numpy.float32), tmp.as_array()
            )
    
        tmp1 = image1.subtract(image2)
        
        numpy.testing.assert_array_almost_equal(
            numpy.zeros(image1.shape, dtype=numpy.float32), tmp1.as_array()
            )
        
        tmp1.fill(2.)
        image1.subtract(image2, out=tmp1)
        
        numpy.testing.assert_array_almost_equal(
            numpy.zeros(image1.shape, dtype=numpy.float32), tmp1.as_array()
            )



class TestGadgetronAcquisitionDataAlgebra(unittest.TestCase):

    def setUp(self):
        os.chdir(examples_data_path('MR'))
        
        self.cwd = os.getcwd()

    
    def tearDown(self):
        #shutil.rmtree(self.cwd)
        pass


    def test_divide_scalar(self):
        os.chdir(self.cwd)
        image1 = mr.AcquisitionData('simulated_MR_2D_cartesian.h5')
        image2 = mr.AcquisitionData('simulated_MR_2D_cartesian.h5')
        image1.sort()
        image2.sort()
        image1.fill(1.)
        image2.fill(2.)
        
        tmp = image1/1.
        numpy.testing.assert_array_equal(image1.as_array(), tmp.as_array())
    
        tmp1 = image1.divide(1.)
        numpy.testing.assert_array_equal(tmp.as_array(), tmp1.as_array())
        
        image1.divide(1., out=image2)
        numpy.testing.assert_array_equal(tmp.as_array(), image2.as_array())
        
    def test_divide_datacontainer(self):
        os.chdir(self.cwd)
        
        image1 = mr.AcquisitionData('simulated_MR_2D_cartesian.h5') 
        image2 = mr.AcquisitionData('simulated_MR_2D_cartesian.h5') 
        image1.sort()
        image2.sort()

        # add 1 because the data contains zeros and divide is not going to be happy
        image1+=1.
        image2+=1.
        tmp = image1/image2

        numpy.testing.assert_array_almost_equal(
            numpy.ones(image1.shape, dtype=numpy.float32), tmp.as_array()
            )
    
        tmp1 = image1.divide(image2)
        numpy.testing.assert_array_almost_equal(
            numpy.ones(image1.shape, dtype=numpy.float32), tmp1.as_array()
            )
        
        tmp1.fill(2.)
        image1.divide(image2, out=tmp1)
        
        numpy.testing.assert_array_almost_equal(
            numpy.ones(image1.shape, dtype=numpy.float32), tmp1.as_array()
            )
        

    def test_multiply_scalar(self):
        os.chdir(self.cwd)
        image1 = mr.AcquisitionData('simulated_MR_2D_cartesian.h5')
        image2 = mr.AcquisitionData('simulated_MR_2D_cartesian.h5')
        image1.sort()
        image2.sort()

        image2.fill(2.)
        
        tmp = image1 * 1.
        numpy.testing.assert_array_equal(image1.as_array(), tmp.as_array())
    
        tmp1 = image1.multiply(1.)
        numpy.testing.assert_array_equal(tmp.as_array(), tmp1.as_array())
        
        image1.multiply(1., out=image2)
        numpy.testing.assert_array_equal(tmp.as_array(), image2.as_array())
        
    def test_multiply_datacontainer(self):
        os.chdir(self.cwd)
        image1 = mr.AcquisitionData('simulated_MR_2D_cartesian.h5')
        image2 = mr.AcquisitionData('simulated_MR_2D_cartesian.h5')
        image1.sort()
        image2.sort()

        image2.fill(1.)
        tmp = image1 * image2

        numpy.testing.assert_array_almost_equal(
            image1.as_array(), tmp.as_array()
            )
    
        tmp1 = image1.multiply(image2)
        numpy.testing.assert_array_almost_equal(
            image1.as_array(), tmp1.as_array()
            )
        
        tmp1.fill(2.)
        image1.multiply(image2, out=tmp1)
        
        numpy.testing.assert_array_almost_equal(
            image1.as_array(), tmp1.as_array()
            )
        
    def test_add_scalar(self):
        os.chdir(self.cwd)
        image1 = mr.AcquisitionData('simulated_MR_2D_cartesian.h5')
        image2 = mr.AcquisitionData('simulated_MR_2D_cartesian.h5')
        image1.sort()
        image2.sort()

        image1.fill(0)
        image2.fill(1)
        
        tmp = image1 + 1.
        numpy.testing.assert_array_equal(image2.as_array(), tmp.as_array())
    
        tmp1 = image1.add(1.)
        numpy.testing.assert_array_equal(tmp.as_array(), tmp1.as_array())
        
        tmp1.fill(0)
        image1.add(1., out=tmp1)
        numpy.testing.assert_array_equal(tmp1.as_array(), image2.as_array())
    
    def test_add_datacontainer(self):
        os.chdir(self.cwd)
        image1 = mr.AcquisitionData('simulated_MR_2D_cartesian.h5')
        image2 = mr.AcquisitionData('simulated_MR_2D_cartesian.h5')
        image1.sort()
        image2.sort()

        image1.fill(0.)
        image2.fill(1.)
        tmp = image1 + image2

        numpy.testing.assert_array_almost_equal(
            numpy.ones(image1.shape, dtype=numpy.float32), tmp.as_array()
            )
    
        tmp1 = image1.add(image2)
        
        numpy.testing.assert_array_almost_equal(
            numpy.ones(image1.shape, dtype=numpy.float32), tmp1.as_array()
            )
        
        tmp1.fill(2.)
        image1.add(image2, out=tmp1)
        
        numpy.testing.assert_array_almost_equal(
            numpy.ones(image1.shape, dtype=numpy.float32), tmp1.as_array()
            )
        
    
    def test_subtract_scalar(self):
        os.chdir(self.cwd)
        image1 = mr.AcquisitionData('simulated_MR_2D_cartesian.h5')
        image2 = mr.AcquisitionData('simulated_MR_2D_cartesian.h5')
        image1.sort()
        image2.sort()

        image1.fill(2)
        image2.fill(1)
        
        tmp = image1 - 1.
        numpy.testing.assert_array_equal(image2.as_array(), tmp.as_array())
    
        tmp1 = image1.subtract(1.)
        numpy.testing.assert_array_equal(tmp.as_array(), tmp1.as_array())
        
        tmp1.fill(0)
        image1.subtract(1., out=tmp1)
        numpy.testing.assert_array_equal(tmp1.as_array(), image2.as_array())
    
    def test_subtract_datacontainer(self):
        os.chdir(self.cwd)
        image1 = mr.AcquisitionData('simulated_MR_2D_cartesian.h5')
        image2 = mr.AcquisitionData('simulated_MR_2D_cartesian.h5')
        image1.sort()
        image2.sort()

        tmp = image1 - image2

        numpy.testing.assert_array_almost_equal(
            numpy.zeros(image1.shape, dtype=numpy.float32), tmp.as_array()
            )
    
        tmp1 = image1.subtract(image2)
        
        numpy.testing.assert_array_almost_equal(
            numpy.zeros(image1.shape, dtype=numpy.float32), tmp1.as_array()
            )
        
        tmp1.fill(2.)
        image1.subtract(image2, out=tmp1)
        
        numpy.testing.assert_array_almost_equal(
            numpy.zeros(image1.shape, dtype=numpy.float32), tmp1.as_array()
            )