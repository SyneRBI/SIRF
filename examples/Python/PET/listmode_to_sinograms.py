'''Listmode-to-sinograms conversion demo.

Usage:
  listmode_to_sinograms [--help | options]

Options:
  -p <path>, --path=<path>     path to data files, defaults to data/examples/PET/mMR
                               subfolder of SIRF root folder
  -l <list>, --list=<list>     listmode file [default: list.l.hdr]
  -o <sino>, --sino=<sino>     output file prefix [default: sinograms]
  -t <tmpl>, --tmpl=<tmpl>     raw data template [default: mMR_template_span11_small.hs]
  -i <int>, --interval=<int>   scanning time interval to convert as string '(a,b)'
                               [default: (0,10)]
  -e <engn>, --engine=<engn>   reconstruction engine [default: STIR]
  -s <stsc>, --storage=<stsc>  acquisition data storage scheme [default: memory]
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

__version__ = '1.0.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

from ast import literal_eval

from sirf.Utilities import show_2D_array

# import engine module
exec('from sirf.' + args['--engine'] + ' import *')


# process command-line options
data_path = args['--path']
if data_path is None:
    # default to data/examples/PET/mMR
    # Note: seem to need / even on Windows
    #data_path = os.path.join(examples_data_path('PET'), 'mMR')
    data_path = examples_data_path('PET') + '/mMR'
prefix = data_path + '/'
list_file = args['--list']
sino_file = args['--sino']
tmpl_file = args['--tmpl']
list_file = existing_filepath(data_path, list_file)
tmpl_file = existing_filepath(data_path, tmpl_file)
interval = literal_eval(args['--interval'])
storage = args['--storage']
show_plot = not args['--non-interactive']


def main():

    # select acquisition data storage scheme
    AcquisitionData.set_storage_scheme(storage)

    # read acquisition data template
    acq_data_template = AcquisitionData(tmpl_file)

    # create listmode-to-sinograms converter object
    lm2sino = ListmodeToSinograms()

    # set input, output and template files
    lm2sino.set_input(list_file)
    lm2sino.set_output_prefix(sino_file)
    # the template is used to specify the sizes of the output sinogram.
    # see the acquisition_data_from_scanner_info demo for an example how to 
    # make your own template file
    lm2sino.set_template(acq_data_template)
    # old way (now just an alternative option)
    # lm2sino.set_template(tmpl_file)

    # set interval
    lm2sino.set_time_interval(interval[0], interval[1])

    # set some flags as examples (the following values are the defaults)
    lm2sino.flag_on('store_prompts')
    lm2sino.flag_off('interactive')

    # set up the converter
    lm2sino.set_up()

    # convert
    lm2sino.process()

    # get access to the sinograms
    acq_data = lm2sino.get_output()
    # copy the acquisition data into a Python array
    acq_array = acq_data.as_array()
    acq_dim = acq_array.shape
    print('acquisition data dimensions: %dx%dx%dx%d' % acq_dim)
    z = acq_dim[1]//2
    if show_plot:
        show_2D_array('Acquisition data', acq_array[0,z,:,:])

    # compute randoms
    print('estimating randoms, please wait...')
    randoms = lm2sino.estimate_randoms()
    rnd_array = randoms.as_array()
    if show_plot:
        show_2D_array('Randoms', rnd_array[0,z,:,:])


try:
    main()
    print('\n=== done with %s' % __file__)

except error as err:
    print('%s' % err.value)
