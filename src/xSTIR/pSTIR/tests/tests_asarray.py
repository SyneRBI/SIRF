'''sirf.STIR.ImageData.asarray tests.
v{version}

Usage:
  tests_asarray.py [--help | options]

Options:
  -f <file>, --file=<file>     raw data file [default: my_forward_projection.hs]
  -p <path>, --path=<path>     path to data files, defaults to data/examples/PET
                               subfolder of SIRF root folder
  -e <engn>, --engine=<engn>   reconstruction engine [default: STIR]
  -r, --record   not used
  -v, --verbose  report each test status

{author}

{licence}
'''
from sirf.STIR import *
from sirf.Utilities import runner, RE_PYEXT, __license__

__version__ = "0.2.3"
__author__ = "Evgueni Ovtchinnikov"

## SyneRBI Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2025 Rutherford Appleton Laboratory STFC
## Copyright 2015 - 2025 University College London.
##
## This is software developed for the Collaborative Computational
## Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
## (http://www.ccpsynerbi.ac.uk/).
##
## Licensed under the Apache License, Version 2.0 (the "License");
##   you may not use this file except in compliance with the License.
##   You may obtain a copy of the License at
##       http://www.apache.org/licenses/LICENSE-2.0
##   Unless required by applicable law or agreed to in writing, software
##   distributed under the License is distributed on an "AS IS" BASIS,
##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
##   See the License for the specific language governing permissions and
##   limitations under the License.

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

#import math
import numpy

from sirf.Utilities import error, examples_data_path, existing_filepath

# import engine module
import importlib
engine = args['--engine']
pet = importlib.import_module('sirf.' + engine)

# process command-line options
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = examples_data_path('PET')
pet.AcquisitionData.set_storage_scheme('memory')


def test_img_asarray(img_data, test):
    try:
        img_asarray = img_data.asarray(copy=False) # view
        print('img_data.asarray() ok')
    except Exception:
        print('data not contiguous, working with its contiguous copy instead...')
        img_data = img_data + 0
        img_asarray = img_data.asarray(copy=False)
    img_as_array = img_data.as_array() # deepcopy to compare with
    diff = img_asarray - img_as_array # must be 0
    #print('norm of img_data.asarray() - img_data.as_array(): %f' % numpy.linalg.norm(diff))
    test.check_if_equal(numpy.linalg.norm(diff), 0)

    img_asarray += 1 # img_data changed too
    img_as_array = img_data.as_array() # must be 0
    diff = img_asarray - img_as_array
    #print('norm of img_data.asarray() - img_data.as_array(): %f' % numpy.linalg.norm(diff))
    test.check_if_equal(numpy.linalg.norm(diff), 0)


def test_main(rec=False, verb=False, throw=True):

    datafile = RE_PYEXT.sub(".txt", __file__)
    test = pTest(datafile, rec, throw=throw)
    test.verbose = verb
    data_path = examples_data_path('PET')

    engine_version = pet.get_engine_version_string()
    print('Using %s version %s as the reconstruction engine' % (engine, engine_version))

    # direct all engine's messages to files
    _ = pet.MessageRedirector('info.txt', 'warn.txt', 'errr.txt')

    # PET acquisition data to be read from this file
    raw_data_file = existing_filepath(data_path, data_file)
    print('raw data: %s' % raw_data_file)
    acq_data = pet.AcquisitionData(raw_data_file)

    # copy the acquisition data into a Python array
    dim = acq_data.dimensions()
    print('data dimensions: %d x %d x %d x %d' % dim)
    acq_array = acq_data.as_array()

    print('\n== Testing asarray() method of pet.AcquisitionData...')

#    acq_asarray = numpy.asarray(acq_data)
    acq_asarray = acq_data.asarray() # same as above
    diff = acq_array - acq_asarray
    #print('norm of acq_data.as_array() - acq_data.asarray(): %f' % numpy.linalg.norm(diff))
    test.check_if_equal(numpy.linalg.norm(diff), 0)

    print('\n== Testing asarray() method of pet.ImageData...')

    print('\n-- testing image from file:')
    img_data = pet.ImageData(existing_filepath(data_path, 'mMR/mu_map.hv')) # contiguous
    test_img_asarray(img_data, test)

    print('\n-- testing image returned by AcquisitionData.create_uniform_image():')
    img_data = acq_data.create_uniform_image(5)
    test_img_asarray(img_data, test)

    print('\n-- testing image constructed from acquisition data:')
    img_data = pet.ImageData(acq_data)
    test_img_asarray(img_data, test)

    return test.failed, test.ntest

if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__)

