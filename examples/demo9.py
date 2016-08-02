'''
GRAPPA reconstruction with the steepest descent step: illustrates
the use of Acquisition Model projections
'''

import argparse
import os
import pylab
import sys

BUILD_PATH = os.environ.get('BUILD_PATH') + '/xGadgetron'
SRC_PATH = os.environ.get('SRC_PATH') + '/xGadgetron/pGadgetron'

sys.path.append(BUILD_PATH)
sys.path.append(SRC_PATH)

from pGadgetron import *

parser = argparse.ArgumentParser(description = \
'''
GRAPPA reconstruction with the steepest descent step: illustrates
the use of Acquisition Model projections
''')
parser.add_argument\
('filename', nargs='?', default = 'testdata.h5', \
 help = 'raw data file name (default: testdata.h5)')
args = parser.parse_args()                                 

def main():

    # acquisitions will be read from an HDF file args.filename
    input_data = MR_Acquisitions(args.filename)

    # pre-process acquisitions
    prep_gadgets = ['NoiseAdjustGadget', 'AsymmetricEchoGadget', \
         'RemoveROOversamplingGadget']
    print('---\n pre-processing acquisitions...')
    preprocessed_data = input_data.process(prep_gadgets)
    pp_norm = preprocessed_data.norm()

    # compute coil sensitivity maps
    csms = MR_CoilSensitivityMaps()
    print('---\n sorting acquisitions...')
    preprocessed_data.sort()
    print('---\n computing sensitivity maps...')
    csms.calculate(preprocessed_data)

    # perform reconstruction
    recon = MR_BasicGRAPPAReconstruction()
    recon.set_input(preprocessed_data)
    print('---\n reconstructing...')
    recon.process()
    # for undersampled acquisition data GRAPPA computes Gfactor images
    # in addition to reconstructed ones
    complex_images, complex_gfactors = recon.get_output()

    # create acquisition model based on the acquisition parameters
    # stored in preprocessed_data and image parameters stored in complex_images
    am = MR_AcquisitionModel(preprocessed_data, complex_images)
    am.set_coil_sensitivity_maps(csms)

    # use the acquisition model (forward projection) to simulate acquisitions
    fwd_data = am.forward(complex_images)
    fwd_norm = fwd_data.norm()

    # compute the difference between real and simulated acquisitions
    diff = fwd_data - preprocessed_data * (fwd_norm/pp_norm)
    rr = diff.norm()/fwd_norm
    print('---\n reconstruction residual norm (rel): %e' % rr)

    # try to improve the reconstruction by the steepest descent step
    g = am.backward(diff)
    w = am.forward(g)
    alpha = (g*g)/(w*w)
    r_complex_imgs = complex_images - g*alpha

    # get real-valued reconstructed and refined images
    print('---\n processing images...')
    images = complex_images.real()
    r_imgs = r_complex_imgs.real()
    nz = images.number()
    print('%d images reconstructed.' % nz)

    # plot images and gfactors
    print('Enter the number of the slice to view it')
    print('(a value outside the range [1 : %d] will stop this loop)'% nz)
    while True:
        s = str(input('slice: '))
        if len(s) < 1:
            break
        z = int(s)
        if z < 1 or z > nz:
            break
        data = images.image_as_array(z - 1)
        rdata = r_imgs.image_as_array(z - 1)
        pylab.figure(z)
        pylab.title('image')
        pylab.imshow(data[0,0,:,:])
        print('Close Figure %d window to continue...' % z)
        pylab.figure(z + nz)
        pylab.title('refined image')
        pylab.imshow(rdata[0,0,:,:])
        print('Close Figure %d window to continue...' % (z + nz))
        pylab.show()

try:
    main()
    print('done')

except error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
