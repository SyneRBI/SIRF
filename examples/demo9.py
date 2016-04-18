import math
import pylab
import sys
import time

sys.path.append('../../build/xGadgetron')
sys.path.append('../pGadgetron')

from pGadgetron import *

try:
    # acquisitions will be read from this HDF file
    input_data = MR_Acquisitions('testdata.h5')

    processed_data = MR_remove_x_oversampling(input_data)

    # perform reconstruction
    br = MR_BasicReconstruction()
    br.set_input(processed_data)
    br.process()
    complex_images = br.get_output()

    # calculate coil sensitivity maps
    csms = MR_CoilSensitivityMaps()
    print('---\n sorting acquisitions...')
    processed_data.sort()
    print('---\n calculating sensitivity maps...')
    csms.calculate(processed_data)

    # create acquisition model based on the acquisition parameters
    # stored in input_data and image parameters stored in complex_images
    am = MR_AcquisitionModel(processed_data, complex_images)

    am.set_coil_sensitivity_maps(csms)

    # post-process reconstructed images
    print('---\n processing images...')
    images = MR_extract_real_images(complex_images)

    # plot obtained images
    for i in range(images.number()):
        data = images.image_as_array(i)
        pylab.figure(i + 1)
        pylab.imshow(data[0,0,:,:])
        pylab.show()

except error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
