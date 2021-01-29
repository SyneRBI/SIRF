#========================================================================
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
import sirf.Reg as reg
from sirf.Utilities import  examples_data_path, TestDataContainerAlgebra, error


class TestNiftiImageDataAlgebra(unittest.TestCase, TestDataContainerAlgebra):

    def setUp(self):
        os.chdir(examples_data_path('Registration'))
        image1 = reg.ImageData('test2.nii.gz')
        image2 = reg.ImageData('test2.nii.gz')
        
        self.image1 = image1
        self.image2 = image2
        self.cwd = os.getcwd()

    
    def tearDown(self):
        #shutil.rmtree(self.cwd)
        pass

    def test_division_by_scalar_zero(self):
        try:
            self.image1 / 0.
            self.assertFalse(True)
        except ZeroDivisionError as derror:
            self.assertTrue(True)
            print ("Caught error", derror)
        except error as derror:
            self.assertTrue(True)
            print ("Caught error", derror)
        
    
    def test_division_by_datacontainer_zero(self):
        try:
            self.image2 *= 0
            tmp = self.image1 / self.image2
            self.assertFalse(True)
        except ZeroDivisionError as derror:
            self.assertTrue(True)
            print ("Caught error", derror)
        except error as derror:
            self.assertTrue(True)
            print ("Caught error", derror)
