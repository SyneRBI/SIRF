import numpy
import os
import shutil
import unittest
import sirf.STIR as pet
from sirf.Utilities import  examples_data_path


class TestSIRFDataContainerAlgebra(unittest.TestCase):
    
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
        # add 1 because the data contains zeros and divide is not going to be happy
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
        # add 1 because the data contains zeros and divide is not going to be happy
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
        # add 1 because the data contains zeros and divide is not going to be happy
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