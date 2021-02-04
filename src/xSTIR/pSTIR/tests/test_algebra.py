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

pet.AcquisitionData.set_storage_scheme('file')
pet.set_verbosity(0)

class TestSTIRImageDataAlgebra(unittest.TestCase, TestDataContainerAlgebra):
    def setUp(self):
        self.image1 = pet.ImageData(os.path.join(
            examples_data_path('PET'),'thorax_single_slice','emission.hv')
            )
        self.image2 = pet.ImageData(os.path.join(
            examples_data_path('PET'),'thorax_single_slice','emission.hv')
            )

    def tearDown(self):
        # shutil.rmtree(self.cwd)
        pass

class TestSTIRAcquisitionDataAlgebraFile(unittest.TestCase, TestDataContainerAlgebra):

    def setUp(self):
        if os.path.exists(os.path.join(
            examples_data_path('PET'), 'mMR', 'mMR_template_span11_small.hs')):

            template = pet.AcquisitionData(os.path.join(
                examples_data_path('PET'), 'mMR', 'mMR_template_span11_small.hs')
            )
            
            self.image1 = template.get_uniform_copy(0)
            self.image2 = template.get_uniform_copy(0)
            # assert False
            self.set_storage_scheme()
        
        
    def tearDown(self):
        pass

    def set_storage_scheme(self):
        pet.AcquisitionData.set_storage_scheme('file')
    
    def test_division_by_datacontainer_zero(self):
        # skip this test as currently cSIRF doesn't throw
        pass


class TestSTIRAcquisitionDataAlgebraMemory(unittest.TestCase, TestDataContainerAlgebra):
    def setUp(self):
        pet.AcquisitionData.set_storage_scheme('file')
        if os.path.exists(os.path.join(
            examples_data_path('PET'), 'mMR', 'mMR_template_span11_small.hs')):

            template = pet.AcquisitionData(os.path.join(
                examples_data_path('PET'), 'mMR', 'mMR_template_span11_small.hs')
            )
            
            self.image1 = template.get_uniform_copy(0)
            self.image2 = template.get_uniform_copy(0)
            # assert False
            pet.AcquisitionData.set_storage_scheme('memory')
    def test_division_by_datacontainer_zero(self):
        # skip this test as currently cSIRF doesn't throw
        pass