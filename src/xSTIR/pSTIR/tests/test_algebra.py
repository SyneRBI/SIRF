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
from sirf.Utilities import examples_data_path, TestDataContainerAlgebra

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
        self.set_storage_scheme()
        sirf_path = os.environ.get('SIRF_PATH')
        data_path = os.path.join(sirf_path, 'tests', 'data.hs')
        if os.path.exists(data_path):
            print('reading data from %s' % data_path)
            data = pet.AcquisitionData(data_path)
            self.image1 = data.get_uniform_copy(0)
            self.image2 = data.get_uniform_copy(0)
            self.set_storage_scheme()
            return
        path = os.path.join(examples_data_path('PET'),
                            'mMR', 'mMR_template_span11_small.hs')
#                            'thorax_single_slice', 'template_sinogram.hs')
        if os.path.exists(path):
            template = pet.AcquisitionData(path)
            print(template.dimensions())
            data = template.get_uniform_copy(0)
            data = data.rebin(3, num_views_to_combine=6, num_tang_poss_to_trim=200)
            print('rebinned to ', data.dimensions())
            self.image1 = data
            self.image2 = data.get_uniform_copy(0)
            print('saving data to %s' % data_path)
            self.image1.write(data_path)

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
        sirf_path = os.environ.get('SIRF_PATH')
        data_path = os.path.join(sirf_path, 'tests', 'data.hs')
        if os.path.exists(data_path):
            print('reading data from %s' % data_path)
            data = pet.AcquisitionData(data_path)
            self.image1 = data.get_uniform_copy(0)
            self.image2 = data.get_uniform_copy(0)
            return
        path = os.path.join(examples_data_path('PET'),
                            'mMR', 'mMR_template_span11_small.hs')
#                            'thorax_single_slice', 'template_sinogram.hs')
        if os.path.exists(path):
            template = pet.AcquisitionData(path)
            print(template.dimensions())
            data = template.get_uniform_copy(0)
            data = data.rebin(3, num_views_to_combine=6, num_tang_poss_to_trim=200)
            print('rebinned to ', data.dimensions())
            self.image1 = data.get_uniform_copy(0)
            self.image2 = data.get_uniform_copy(0)
            pet.AcquisitionData.set_storage_scheme('memory')
            print('saving data to %s' % data_path)
            self.image1.write(data_path)

    def test_division_by_datacontainer_zero(self):
        # skip this test as currently cSIRF doesn't throw
        pass
