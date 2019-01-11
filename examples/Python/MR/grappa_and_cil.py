'''
Example of reconstruction using CCPi/CIL FISTA algorithm 

Usage:
  grappa_and_cil.py [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: simulated_MR_2D_cartesian_Grappa2.h5]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/MR
                              subfolder of SIRF root folder
'''

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC.
## Copyright 2015 - 2017 University College London.
##
## This is software developed for the Collaborative Computational
## Project in Positron Emission Tomography and Magnetic Resonance imaging
## (http://www.ccppetmr.ac.uk/).
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

from sirf.Utilities import existing_filepath
from sirf.Utilities import error
from sirf.Utilities import show_3D_array

from sirf.Gadgetron import petmr_data_path
from sirf.Gadgetron import AcquisitionData
from sirf.Gadgetron import AcquisitionModel
from sirf.Gadgetron import AcquisitionDataProcessor
from sirf.Gadgetron import CartesianGRAPPAReconstructor
from sirf.Gadgetron import CoilSensitivityData

from ccpi.optimisation.funcs import Norm2sq
from ccpi.optimisation.funcs import ZeroFun
from ccpi.optimisation.algs import FISTA
from ccpi.optimisation.ops import PowerMethodNonsquare

import numpy


# process command-line options
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = petmr_data_path('mr')

def main():

    # locate the input data file
    input_file = existing_filepath(data_path, data_file)

    # acquisition data will be read from an HDF file input_data
    acq_data = AcquisitionData(input_file)

    # pre-process acquisition data
    print('---\n pre-processing acquisition data...')
    # the following function runs a GadgetChain data processor
    # preprocessed_data = preprocess_acquisition_data(acq_data)
    gadget_chain = ['NoiseAdjustGadget', \
                    'AsymmetricEchoAdjustROGadget', \
                    'RemoveROOversamplingGadget']
    ap = AcquisitionDataProcessor( gadget_chain )
    ap.set_input( acq_data )
    ap.process()
    preprocessed_data = ap.get_output()


    # perform reconstruction
    recon = CartesianGRAPPAReconstructor()
    recon.set_input( preprocessed_data )
    recon.compute_gfactors( False )
    print('---\n reconstructing...')
    recon.process()
    # for undersampled acquisition data GRAPPA computes Gfactor images
    # in addition to reconstructed ones
    image_data = recon.get_output()

    # compute coil sensitivity maps
    csms = CoilSensitivityData()
    print('---\n sorting acquisition data...')
    preprocessed_data.sort()
    print('---\n computing sensitivity maps...')
    csms.calculate( preprocessed_data )

    # create acquisition model based on the acquisition parameters
    # stored in preprocessed_data and image parameters stored in
    # image_data
    acq_model = AcquisitionModel( preprocessed_data, image_data )
    acq_model.set_coil_sensitivity_maps( csms )

    # use the acquisition model (forward projection) to simulate acquisition data
    simulated_data = acq_model.forward( image_data )

    # USE FISTA
    little_value = 1e-2
    x_init = image_data.copy() * little_value
    # x_init.fill(numpy.random.randn(*image_data.as_array().shape))
    # x_init.fill(numpy.zeros(numpy.shape(image_data.as_array().shape))+little_value)
    norm2sq = Norm2sq( A = acq_model , b = simulated_data , c = 1)
    no_regulariser = ZeroFun()
    options = {'tol': 1e-4, 'iter': 10, 'memopt':False}

    # calculate Lipschitz constant
    # x_init.fill(numpy.random.randn(*x_init.as_array().shape))

    norm2sq.L = PowerMethodNonsquare( acq_model , numiters = 1 , x0 = x_init) [0]
    #norm2sq.L = 0.5

    print ("Lipschitz " , norm2sq.L)
    x_fista0, it0, timing0, criter0 = FISTA(x_init, norm2sq, no_regulariser ,  opt=options)

    image_array = image_data.as_array()
    refined_image_array = x_fista0.as_array()

    # show reconstructed and refined images
    title = 'Reconstructed image data (magnitude)'
    show_3D_array(abs(image_array), suptitle = title, label = 'slice', \
                  xlabel = 'samples', ylabel = 'readouts', show = False)
    title = 'FISTA image data (magnitude)'
    show_3D_array(abs(refined_image_array), suptitle = title, label = 'slice', \
                  xlabel = 'samples', ylabel = 'readouts')

try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)
