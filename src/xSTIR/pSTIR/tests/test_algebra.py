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

import os
import unittest
import sirf.STIR as pet
from sirf.Utilities import examples_data_path, DataContainerAlgebraTests

pet.AcquisitionData.set_storage_scheme('file')
pet.set_verbosity(0)

class TestSTIRImageDataAlgebra(unittest.TestCase, DataContainerAlgebraTests):
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

class TestSTIRAcquisitionDataAlgebraFile(unittest.TestCase, DataContainerAlgebraTests):

    def setUp(self):
        path = os.path.join(
            examples_data_path('PET'), 'thorax_single_slice', 'template_sinogram.hs')
        if os.path.exists(path):
            template = pet.AcquisitionData(path)
            self.image1 = template.get_uniform_copy(0)
            self.image2 = template.get_uniform_copy(0)
            self.set_storage_scheme()
        
        
    def tearDown(self):
        pass

    def set_storage_scheme(self):
        pet.AcquisitionData.set_storage_scheme('file')
    
    def test_division_by_datacontainer_zero(self):
        # skip this test as currently cSIRF doesn't throw
        pass


class TestSTIRAcquisitionDataAlgebraMemory(unittest.TestCase, DataContainerAlgebraTests):
    def setUp(self):
        pet.AcquisitionData.set_storage_scheme('file')
        path = os.path.join(
            examples_data_path('PET'), 'thorax_single_slice', 'template_sinogram.hs')
        if os.path.exists(path):
            template = pet.AcquisitionData(path)
            self.image1 = template.get_uniform_copy(0)
            self.image2 = template.get_uniform_copy(0)
            pet.AcquisitionData.set_storage_scheme('memory')
    def test_division_by_datacontainer_zero(self):
        # skip this test as currently cSIRF doesn't throw
        pass

class TestSTIRAcquisitionDataSubsetAlgebraMemory(unittest.TestCase, DataContainerAlgebraTests):
    def setUp(self):
        pet.AcquisitionData.set_storage_scheme('memory')
        path = os.path.join(
            examples_data_path('PET'), 'thorax_single_slice', 'template_sinogram.hs')
        if os.path.exists(path):
            template = pet.AcquisitionData(path)
            prompts = template.get_uniform_copy(0)
            # create a staggered list of views for 2 subsets
            views = prompts.dimensions()[2]
            indices = list(range(views))
            num_batches = 2
            batches = [indices[i::num_batches] for i in range(num_batches)]
            
            self.image1 = prompts.get_subset(batches[0])
            self.image2 = prompts.get_subset(batches[0])
            
    def test_division_by_datacontainer_zero(self):
        # skip this test as currently cSIRF doesn't throw
        pass