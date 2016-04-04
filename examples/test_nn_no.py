import math
import pylab
import sys
import time

sys.path.append('../../build/xGadgetron')
sys.path.append('../pGadgetron')

from pGadgetron import *

try:
    # acquisitions will be read from this HDF file
    input_data = ISMRMRDAcquisitions('nn_no.h5')

    # perform reconstruction
    recon = SimpleReconstructionProcessor()
    recon.set_input(input_data)
    print('reconstructing...')
    recon.process()
    interim_images = recon.get_output()

    csms = MRCoilSensitivityMaps()

    print('ordering acquisitions...')
    input_data.order()

    print('computing sensitivity maps...')
    csms.compute(input_data)

    # create acquisition model based on the acquisition parameters
    # stored in input_data and image parameters stored in interim_images
    am = AcquisitionModel(input_data, interim_images)

    am.set_coil_sensitivity_maps(csms)

    # use the acquisition model (forward projection) to produce acquisitions
    acqs = am.forward(interim_images)

    print('---\n their forward projection norm %e' % acqs.norm())

    # compute the difference between real and modelled acquisitions:
    #   diff = acqs - P acqs,
    # where P is the orthogonal projector onto input_data
    a = -acqs.dot(input_data) / input_data.dot(input_data)
    b = 1.0
    diff = AcquisitionsContainer.axpby(a, input_data, b, acqs)
    rr = diff.norm()/acqs.norm()
    print('---\n reconstruction residual norm (rel): %e' % rr)

    # post-process reconstructed images
    print('processing images...')
    images = MR_extract_real_images(interim_images)

    # plot obtained images
    for i in range(images.number()):
        data = images.image_as_array(i)
        pylab.figure(i + 1)
        pylab.imshow(data[0,0,:,:])
        pylab.show()

except error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
