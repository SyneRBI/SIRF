'''
Upper-level demo, GRAPPA reconstruction of undersampled data.

Usage:
  undersampled_recon.py [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: simulated_MR_2D_cartesian_Grappa2.h5]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/MR
                              subfolder of SIRF root folder
  -e <engn>, --engine=<engn>  reconstruction engine [default: Gadgetron]
'''

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
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

try:
    import pylab
    HAVE_PYLAB = True
except:
    HAVE_PYLAB = False

# import engine module
exec('from p' + args['--engine'] + ' import *')

def main():

    # locate the input data file
    data_path = args['--path']
    if data_path is None:
        data_path = mr_data_path()
    input_file = existing_filepath(data_path, args['--file'])

    # acquisitions will be read from an HDF file input_file
    input_data = AcquisitionData(input_file)
    if not input_data.is_undersampled():
        print('this demo needs undersampled raw data')
        return

    # pre-process acquisitions
    print('---\n pre-processing acquisitions...')
    preprocessed_data = preprocess_acquisitions(input_data)

    # perform reconstruction
    recon = CartesianGRAPPAReconstruction()
    # for undersampled acquisition data GRAPPA will compute Gfactor images
    # in addition to reconstructed ones
    recon.compute_gfactors(True)
    recon.set_input(preprocessed_data)
    print('---\n reconstructing...')
    recon.process()
    image = recon.get_output('image')
    gfact = recon.get_output('gfactor')
    idata = abs(image.as_array())
    gdata = abs(gfact.as_array())

    nz = idata.shape[0]
    # plot image and gfactor slices
    while HAVE_PYLAB:
        print('---\n Enter the slice number to view it.')
        print(' A value outside the range [1 : %d] will stop this loop.'% nz)
        s = str(input('slice: '))
        if len(s) < 1:
            break
        z = int(s)
        if z < 1 or z > nz:
            break
        pylab.figure(z)
        pylab.title('image')
        pylab.imshow(idata[z - 1,:,:])
        print('Close Figure %d window to continue...' % z)
        pylab.figure(z + nz)
        pylab.title('G factor')
        pylab.imshow(gdata[z - 1,:,:])
        print('Close Figure %d window to continue...' % (z + nz))
        pylab.show()

try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)
