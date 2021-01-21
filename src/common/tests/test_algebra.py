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

    def test_divide(self):
        os.chdir(self.cwd)
        image1 = pet.ImageData('emission.hv')
        image2 = pet.ImageData('emission.hv')
        image1.fill(1.)
        image2.fill(2.)
        print (image1.shape, image2.shape)
        
        
        tmp = image1/1.
        numpy.testing.assert_array_equal(image1.as_array(), tmp.as_array())
    
        tmp1 = image1.divide(1.)
        numpy.testing.assert_array_equal(tmp.as_array(), tmp1.as_array())
        
        image1.divide(1., out=image2)
        numpy.testing.assert_array_equal(tmp.as_array(), image2.as_array())
        

    def test_multiply(self):
        os.chdir(self.cwd)
        image1 = pet.ImageData('emission.hv')
        image2 = pet.ImageData('emission.hv')
        
        image2.fill(2.)
        print (image1.shape, image2.shape)
        
        tmp = image1 * 1.
        numpy.testing.assert_array_equal(image1.as_array(), tmp.as_array())
    
        tmp1 = image1.multiply(1.)
        numpy.testing.assert_array_equal(tmp.as_array(), tmp1.as_array())
        
        image1.multiply(1., out=image2)
        numpy.testing.assert_array_equal(tmp.as_array(), image2.as_array())
        
    
    def test_add(self):
        os.chdir(self.cwd)
        image1 = pet.ImageData('emission.hv')
        image2 = pet.ImageData('emission.hv')
        image1.fill(0)
        image2.fill(1)
        print (image1.shape, image2.shape)
        
        tmp = image1 + 1.
        numpy.testing.assert_array_equal(image2.as_array(), tmp.as_array())
    
        tmp1 = image1.add(1.)
        numpy.testing.assert_array_equal(tmp.as_array(), tmp1.as_array())
        
        tmp1.fill(0)
        image1.add(1., out=tmp1)
        numpy.testing.assert_array_equal(tmp1.as_array(), image2.as_array())
    
    def test_subtract(self):
        os.chdir(self.cwd)
        image1 = pet.ImageData('emission.hv')
        image2 = pet.ImageData('emission.hv')
        image1.fill(2)
        image2.fill(1)
        print (image1.shape, image2.shape)
        
        tmp = image1 - 1.
        numpy.testing.assert_array_equal(image2.as_array(), tmp.as_array())
    
        tmp1 = image1.subtract(1.)
        numpy.testing.assert_array_equal(tmp.as_array(), tmp1.as_array())
        
        tmp1.fill(0)
        image1.subtract(1., out=tmp1)
        numpy.testing.assert_array_equal(tmp1.as_array(), image2.as_array())