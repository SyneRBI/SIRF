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

import numpy
import os
import shutil
import unittest
import sirf.STIR as pet
from sirf.Utilities import examples_data_path, TestDataContainerAlgebra, error

pet.AcquisitionData.set_storage_scheme('memory')

class TestSTIRImageDataAlgebra(unittest.TestCase, TestDataContainerAlgebra):
    def setUp(self):
        self.set_storage_scheme()
        self.image1 = pet.ImageData(os.path.join(
            examples_data_path('PET'),'thorax_single_slice','emission.hv')
            )
        self.image2 = pet.ImageData(os.path.join(
            examples_data_path('PET'),'thorax_single_slice','emission.hv')
            )
        
        
    def tearDown(self):
        # shutil.rmtree(self.cwd)
        pass

    def set_storage_scheme(self):
        pet.AcquisitionData.set_storage_scheme('memory')

class TestSTIRAcquisitionDataAlgebra(unittest.TestCase, TestDataContainerAlgebra):

    def setUp(self):
        self.set_storage_scheme()
        
        self.image1 = pet.ImageData(os.path.join(
            examples_data_path('PET'),'thorax_single_slice','emission.hv')
            )
        self.image1 = self.image1.get_uniform_copy(0)
        
        self.image2 = pet.ImageData(os.path.join(
            examples_data_path('PET'),'thorax_single_slice','emission.hv')
            )
        self.image2 = self.image2.get_uniform_copy(0)
        
        
    def tearDown(self):
        pass

    def set_storage_scheme(self):
        pet.AcquisitionData.set_storage_scheme('memory')

class TestSTIRImageDataAlgebraFile(TestSTIRImageDataAlgebra):
    def set_storage_scheme(self):
        pet.AcquisitionData.set_storage_scheme('file')

class TestSTIRAcquisitionDataAlgebraFile(TestSTIRAcquisitionDataAlgebra):

    def set_storage_scheme(self):
        pet.AcquisitionData.set_storage_scheme('file')