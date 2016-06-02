'''
Upper level 3-steps GRAPPA reconstruction demo.
'''

import os
import pylab
import sys
import time

BUILD_PATH = os.environ.get('BUILD_PATH') + '/xGadgetron'
SRC_PATH = os.environ.get('SRC_PATH') + '/xGadgetron/pGadgetron'

sys.path.append(BUILD_PATH)
sys.path.append(SRC_PATH)

from pGadgetron import *

try:
    # acquisitions will be read from this HDF file
    file = str(input('raw data file (with apostrophys in Python2.*): '))
    input_data = MR_Acquisitions(file)

    # pre-process acquisitions
    prep_gadgets = ['NoiseAdjustGadget', 'AsymmetricEchoGadget', \
         'RemoveROOversamplingGadget']
    acq_proc = AcquisitionsProcessor(prep_gadgets)
    print('pre-processing acquisitions...')
    preprocessed_data = acq_proc.process(input_data)

    csms = MR_CoilSensitivityMaps()
    print('---\n sorting acquisitions...')
    preprocessed_data.sort()
    print('---\n computing sensitivity maps...')
    csms.calculate(preprocessed_data)

    recon = MR_BasicGRAPPAReconstruction()
    # connect to input data
    recon.set_input(preprocessed_data)
    # perform reconstruction
    print('reconstructing...')
    recon.process()
    # get reconstructed images
    complex_images = recon.get_output()

    # post-process reconstructed images
    print('processing images...')
    images = MR_extract_real_images(complex_images)

    # for undersampled acquisition data GRAPPA computes Gfactor images
    # in addition to reconstructed ones
    nt = images.types() # 1 for fully sampled, 2 for undersampled
    ni = images.number()
    nz = ni/nt
    print('%d images reconstructed.' % nz)

    # plot obtained images
    print('Enter z-coordinate of the slice to view it')
    print('(a value outside the range [0 : %d] will stop this loop)'%(nz - 1))
    while True:
        s = str(input('z-coordinate: '))
        if len(s) < 1:
            break
        z = int(s)
        if z < 0 or z >= nz:
            break
        for j in range(nt):
            i = nt*z + j
            data = images.image_as_array(i)
            pylab.figure(i + 1)
            pylab.imshow(data[0,0,:,:])
            print('Close Figure %d window to continue...' % (i + 1))
        pylab.show()

except error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
