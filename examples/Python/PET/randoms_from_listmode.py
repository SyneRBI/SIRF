'''randoms-from-listmode estimation demo

Usage:
  randoms_from_listmode [--help | options]

This demo illustrates how to get a good estimate of the randoms from a list mode file.
It compares the result with the original delayed coincidences.

This demo will therefore only work with list mode data where the delayeds are stored.

Options:
  -p <path>, --path=<path>     path to data files, defaults to data/examples/PET/mMR
                               subfolder of SIRF root folder
  -l <list>, --list=<list>     listmode file [default: list.l.hdr]
  -o <sino>, --sino=<sino>     sinograms file prefix [default: sinograms]
  -r <rand>, --rand=<rand>     randoms file [default: randoms]
  -t <tmpl>, --tmpl=<tmpl>     raw data template [default: mMR_template_span11_small.hs]
  -i <int>, --interval=<int>   scanning time interval to convert as string '(a,b)'
                               (no space after comma) [default: (0,100)]
  -e <engn>, --engine=<engn>   reconstruction engine [default: STIR]
  -s <stsc>, --storage=<stsc>  acquisition data storage scheme [default: file]
  --non-interactive            do not show plots
'''
## SyneRBI Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2018 - 2019 Rutherford Appleton Laboratory STFC
## Copyright 2018 University College London.
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
import importlib
from ast import literal_eval

import numpy as np
from docopt import docopt
from sirf.Utilities import examples_data_path, existing_filepath, show_3D_array

try:
    import pylab
    HAVE_PYLAB = True
except RuntimeWarning:
    HAVE_PYLAB = False
__version__ = '0.1.0'


def main(argv):
    # process command-line options
    args = docopt(__doc__, version=__version__, argv=argv)
    data_path = args['--path']
    if data_path is None:
        # default to data/examples/PET/mMR
        # Note: seem to need / even on Windows
        #data_path = os.path.join(examples_data_path('PET'), 'mMR')
        data_path = examples_data_path('PET') + '/mMR'
    list_file = args['--list']
    sino_file = args['--sino']
    rand_file = args['--rand']
    tmpl_file = args['--tmpl']
    list_file = existing_filepath(data_path, list_file)
    tmpl_file = existing_filepath(data_path, tmpl_file)
    interval = literal_eval(args['--interval'])
    storage = args['--storage']
    show_plot = not args['--non-interactive'] and HAVE_PYLAB
    # import engine module
    engine = args['--engine']
    pet = importlib.import_module('sirf.' + engine)

    # direct all engine's messages to files
    _ = pet.MessageRedirector('info.txt', 'warn.txt', 'errr.txt')

    # select acquisition data storage scheme
    scheme = pet.AcquisitionData.get_storage_scheme()
    pet.AcquisitionData.set_storage_scheme(storage)

    # create listmode-to-sinograms converter object
    lm2sino = pet.ListmodeToSinograms()

    # set input, output and template files
    lm2sino.set_input(list_file)
    lm2sino.set_output_prefix(sino_file)
    lm2sino.set_template(tmpl_file)

    # set interval
    lm2sino.set_time_interval(interval[0], interval[1])

    # set flags such that we only get the delayed coincidences
    lm2sino.flag_on('store_delayeds')
    lm2sino.flag_off('store_prompts')

    # set up the converter
    lm2sino.set_up()

    # convert
    lm2sino.process()

    # get access to the sinograms
    delayeds_acq_data = lm2sino.get_output()

    # estimate the randoms from the delayeds via Maximum Likelihood estimation
    # This will take at least a few seconds
    randoms_estimate_acq_data = lm2sino.estimate_randoms();
    randoms_estimate_acq_data.write(rand_file)

    # copy the acquisition data into Python arrays
    delayeds_acq_array = delayeds_acq_data.as_array()
    randoms_estimate_acq_array = randoms_estimate_acq_data.as_array()
    acq_dim = delayeds_acq_array.shape
    print('acquisition data dimensions: %dx%dx%dx%d' % acq_dim)
    print('The total number of delayed coincidences and estimated randoms have to be very similar.')
    print('Let us check this:')
    print('total delayeds: %.1f, total estimated randoms: %.1f' % (delayeds_acq_array.sum(), randoms_estimate_acq_array.sum()))
    print('Max values should be somewhat similar, but this depends on statistics of course.')
    print('max delayeds: %f, max estimated randoms: %f' % (delayeds_acq_array.max(), randoms_estimate_acq_array.max()))

    print('A single sinogram (this will look very different for noisy data)')
    z = acq_dim[1]//2
    if show_plot:
        show_3D_array(np.stack((delayeds_acq_array[0,z,:,:], randoms_estimate_acq_array[0,z,:,:])), titles=('raw delayeds', ' estimated randoms'))
        pylab.show()

    if scheme != storage:
        pet.AcquisitionData.set_storage_scheme(scheme)


if __name__ == "__main__":
    main(None)
