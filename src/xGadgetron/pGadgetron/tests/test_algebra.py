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
import sirf.Gadgetron as mr
from sirf.Utilities import  examples_data_path, TestDataContainerAlgebra



class TestGadgetronAcquisitionDataAlgebra(unittest.TestCase, TestDataContainerAlgebra):

    def setUp(self):
        image1 = mr.AcquisitionData(os.path.join
            (examples_data_path('MR'),'simulated_MR_2D_cartesian.h5')
        )
        image2 = mr.AcquisitionData(os.path.join
            (examples_data_path('MR'),'simulated_MR_2D_cartesian.h5')
        )

        image1.sort()
        image2.sort()
        self.image1 = image1
        self.image2 = image2

    def tearDown(self):
        pass




class TestGadgetronImageDataAlgebra(unittest.TestCase, TestDataContainerAlgebra):

    def setUp(self):
        acq_data = mr.AcquisitionData(os.path.join
            (examples_data_path('MR'),'simulated_MR_2D_cartesian.h5')
        )
        # ad1.sort()
        preprocessed_data = mr.preprocess_acquisition_data(acq_data)
        recon = mr.FullySampledReconstructor()
        recon.set_input(preprocessed_data)
        recon.process()
        image_data = recon.get_output()
        
        self.image1 = image_data
        self.image2 = image_data * 1
        
    def tearDown(self):
        #shutil.rmtree(self.cwd)
        pass

