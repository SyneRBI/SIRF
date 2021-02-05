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
import sirf.Reg as reg
from sirf.Utilities import  examples_data_path, TestDataContainerAlgebra


class TestNiftiImageDataAlgebra(unittest.TestCase, TestDataContainerAlgebra):

    def setUp(self):
        image1 = reg.ImageData(os.path.join(
            examples_data_path('Registration'),'test2.nii.gz')
        )
        image2 = reg.ImageData(os.path.join(
            examples_data_path('Registration'),'test2.nii.gz')
        )
        
        self.image1 = image1
        self.image2 = image2

    def tearDown(self):
        #shutil.rmtree(self.cwd)
        pass
