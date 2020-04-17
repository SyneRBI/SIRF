'''
GRAPPA reconstruction with the steepest descent refinement: illustrates
the use of Acquisition Model projections

Usage:
  grappa_and_steepest_descent.py [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: simulated_MR_2D_cartesian_Grappa2.h5]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/MR
                              subfolder of SIRF root folder
  -i <iter>, --iter=<iter>    the number of steepest descent iterations
                              [default: 5]
  -z <zdim>, --zdim=<zdim>    dimension of 3D image data array to be used as z
                              when displaying [default: 0]
  -s <slce>, --slice=<slce>   image slice to display [default: 0]
  -e <engn>, --engine=<engn>  reconstruction engine [default: Gadgetron]
  -o <file>, --output=<file>  images output file
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

# import engine module
exec('from p' + args['--engine'] + ' import *')

# process command-line options
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = examples_data_path('MR')
niter = int(args['--iter'])
slc = int(args['--slice'])
output_file = args['--output']
if slc < 0:
    slc = None
zdim = args['--zdim']
if zdim == 1:
    zyx = (1, 0, 2)
elif zdim == 2:
    zyx = (2, 1, 0)
else:
    zyx = None

def main():

    # locate the input data file
    input_file = existing_filepath(data_path, data_file)

    # acquisition data will be read from an HDF file input_data
    acq_data = AcquisitionData(input_file)

    # pre-process acquisition data
    print('---\n pre-processing acquisition data...')
    preprocessed_data = preprocess_acquisition_data(acq_data)
    preprocessed_data_norm = preprocessed_data.norm()

    # perform reconstruction
    recon = CartesianGRAPPAReconstructor()
    recon.set_input(preprocessed_data)
    recon.compute_gfactors(False)
    print('---\n reconstructing...')
    recon.process()
    # for undersampled acquisition data GRAPPA computes Gfactor images
    # in addition to reconstructed ones
    image_data = recon.get_output()
    title = 'Reconstructed image data (magnitude)'
    image_data.show(zyx=zyx, slice=slc, title=title, \
                    postpone=(niter > 0), cmap=None)
    if niter < 1:
        return

    # compute coil sensitivity maps
    csms = CoilSensitivityData()
    print('---\n sorting acquisition data...')
    preprocessed_data.sort()
    print('---\n computing sensitivity maps...')
    csms.calculate(preprocessed_data)

    # create acquisition model based on the acquisition parameters
    # stored in preprocessed_data and image parameters stored in
    # image_data
    acq_model = AcquisitionModel(preprocessed_data, image_data)
    acq_model.set_coil_sensitivity_maps(csms)

    res = numpy.ndarray((niter,))
    scale = 1 # reconstructed image and phantom may have different scale,
              # and so would forward projection and pre-processed raw data
              # hence a proper scaling of raw data is needed
    for i in range(niter):
        it = i + 1
        simulated_data = acq_model.forward(image_data)
        simulated_data_norm = simulated_data.norm()
        if i == 0:
            scale = preprocessed_data_norm/simulated_data_norm
        residual = simulated_data - preprocessed_data/scale
        res[i] = residual.norm()/preprocessed_data_norm
        print('---- iteration %d, residual norm: %.3e' % (it, res[i]))
        # compute gradient
        grad = acq_model.backward(residual)
        # compute locally optimal steepest descent step
        w = acq_model.forward(grad)
        tau = (grad.dot(grad))/(w.dot(w))
        image_data = image_data - grad * tau
        if i%10 == 0 or i == niter - 1:
            it = i + 1
            title = 'Steepest-descent-refined image data, iteration %d' % it
            image_data.show(zyx=zyx, slice=slc, title=title, cmap=None, \
                            postpone=(i < niter - 1))

    if niter > 1:
        import pylab
        pylab.figure()
        pylab.plot(numpy.arange(1, niter + 1, 1), res)
        pylab.grid()
        pylab.title('residual norm')
        pylab.show()

    if output_file is not None:
      # write images to a new group in args.output
      # named after the current date and time
      time_str = time.asctime()
      print('writing to %s' % output_file)
      image_data.write(output_file) #, time_str)

try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)
