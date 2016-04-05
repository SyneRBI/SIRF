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

    processed_data = MR_remove_x_oversampling(input_data)

    # perform reconstruction
    recon = SimpleReconstructionProcessor()
    recon.set_input(processed_data)
    recon.process()
    complex_images = recon.get_output()

    csms = MRCoilSensitivityMaps()

    # coil sensitivity maps can be read from a file
##    csm_file = str(input('csm file: '))
##    print('reading sensitivity maps...')
##    csms.read(csm_file)
    # or computed
    print('---\n ordering acquisitions...')
    input_data.order()
    print('---\n computing sensitivity maps...')
    csms.compute(input_data)

    # create acquisition model based on the acquisition parameters
    # stored in input_data and image parameters stored in interim_images
    am = AcquisitionModel(input_data, complex_images)

    am.set_coil_sensitivity_maps(csms)

    # use the acquisition model (forward projection) to produce acquisitions
    acqs = am.forward(complex_images)

    # compute the difference between real and modelled acquisitions:
    #   diff = acqs - P acqs,
    # where P is the orthogonal projector onto input_data
    a = (acqs * input_data) / (input_data * input_data)
    diff = acqs - a * input_data
    rr = diff.norm()/acqs.norm()
    print('---\n reconstruction residual norm (rel): %e' % rr)

    # apply the adjoint model (backward projection)
    imgs = am.backward(diff)

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
