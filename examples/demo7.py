'''
Upper-level demo, GRAPPA reconstruction of undersampled data.
'''

import math
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
    print('---\n pre-processing acquisitions...')
    preprocessed_data = acq_proc.process(input_data)

    # perform reconstruction
    recon = MR_BasicGRAPPAReconstruction()
    recon.set_input(preprocessed_data)
    print('---\n reconstructing...')
    recon.process()
    output = recon.get_output()
    # for undersampled acquisition data GRAPPA computes Gfactor images
    # in addition to reconstructed ones
    complex_images = output.select(2)
    complex_gfactors = output.select(2, 1)

    # get real-valued reconstructed images and gfactors
    print('---\n processing images...')
    images = MR_extract_real_images(complex_images)
    gfactors = MR_extract_real_images(complex_gfactors)

    nz = images.number()
    print('%d images reconstructed.' % nz)

    # plot images and gfactors
    print('Enter z-coordinate of the slice to view it')
    print('(a value outside the range [0 : %d] will stop this loop)'%(nz - 1))
    while True:
        s = str(input('z-coordinate: '))
        if len(s) < 1:
            break
        z = int(s)
        if z < 0 or z >= nz:
            break
        i = z
        data = images.image_as_array(i)
        gdata = gfactors.image_as_array(i)
        pylab.figure(i + 1)
        pylab.title('image')
        pylab.imshow(data[0,0,:,:])
        print('Close Figure %d window to continue...' % (i + 1))
        pylab.figure(i + nz + 1)
        pylab.title('G factor')
        pylab.imshow(gdata[0,0,:,:])
        print('Close Figure %d window to continue...' % (i + nz + 1))
        pylab.show()

    print('done')

except error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
