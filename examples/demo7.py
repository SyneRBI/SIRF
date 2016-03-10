import math
import pylab
import sys
import time

sys.path.append('../../build/xGadgetron')
sys.path.append('../pGadgetron')

from pGadgetron import *

try:
    # acquisitions will be read from this HDF file
    input_data = ISMRMRDAcquisitions('testdata.h5')

    # pre-process acquisition data
    print('processing acquisitions...')
    interim_data = MR_remove_x_oversampling(input_data)

    # perform reconstruction
    recon = SimpleReconstructionProcessor()
    recon.set_input(interim_data)
    print('reconstructing...')
    recon.process()
    interim_images = recon.get_output()

    # post-process reconstructed images
    print('processing images...')
    images = MR_extract_real_images(interim_images)

    # plot obtained images
    for i in range(images.number()):
        data = images.image_as_array(i)
        pylab.figure(i + 1)
        pylab.imshow(data[:,:,0,0])
        pylab.show()

except error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
