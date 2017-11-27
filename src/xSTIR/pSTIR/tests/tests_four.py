'''pSTIR test set 4.

Acquisitions and images algebra tests.

CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
2017 Casper da Costa-Luis

This is software developed for the Collaborative Computational
Project in Positron Emission Tomography and Magnetic Resonance imaging
(http://www.ccppetmr.ac.uk/).

Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at
      http://www.apache.org/licenses/LICENSE-2.0
  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
'''
import math
from pSTIR import *

__version__ = '0.2.0'


def test_main(rec=False, verb=False):
    test = pTest('test4.txt', rec)
    test.verbose = verb

    msg_red = MessageRedirector()

    data_path = petmr_data_path('pet')
    raw_data_file = existing_filepath(data_path, 'my_forward_projection.hs')
    acq_data = AcquisitionData(raw_data_file)

    print('Checking acquisition data algebra:')
    new_acq_data = acq_data.clone()
    diff = new_acq_data - acq_data
    test.check(diff.norm())
    test.check(1 - math.sqrt(acq_data*acq_data)/acq_data.norm())
    new_acq_data = acq_data * 10.0
    test.check(1 - 10*acq_data.norm()/new_acq_data.norm())

    print('Checking images algebra:')
    image_data = acq_data.create_uniform_image(10.0)
    diff = image_data.clone() - image_data
    test.check(diff.norm())
    test.check(1 - math.sqrt(image_data*image_data)/image_data.norm())
    new_image_data = image_data*10
    test.check(1 - 10*image_data.norm()/new_image_data.norm())
