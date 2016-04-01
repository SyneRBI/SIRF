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

    print('---\n acquisition data norm: %e' % input_data.norm())

    interim_data = MR_remove_x_oversampling(input_data)

    print('---\n processed acquisition data norm: %e' % interim_data.norm())

    # perform reconstruction
    recon = SimpleReconstructionProcessor()
    recon.set_input(interim_data)
    recon.process()
    interim_images = recon.get_output()

    print('---\n reconstructed images norm: %e' % interim_images.norm())

##    #s = str(input('csm file: '))
##    s = 'csm_testdata.h5'
##    #s = 'csm_opismrmrd.h5'
##    csms = MRCoilSensitivityMaps(s)
    csms = MRCoilSensitivityMaps()

    print('ordering acquisitions...')
    input_data.order()

    print('computing sensitivity maps...')
    csms.compute(input_data)

    nz = csms.number()
    print('%d slices' % nz)

    maxv = 0
    for z in range(nz):
        data = csms.csm_as_array(z)
        minvz = numpy.amin(data)
        maxvz = numpy.amax(data)
        if z == 0:
            minv = minvz
        else:
            minv = min(minvz, minv)
        maxv = max(maxvz, maxv)
    print(minv, maxv)

    while True:
        s = str(input('enter z-coordinate: '))
        if len(s) < 1:
            break
        z = int(s)
        if z < 0 or z >= nz:
            break
        data = csms.csm_as_array(z)/maxv
        shape = data.shape
        nc = shape[0]
        for i in range(nc):
            pylab.figure(z*nc + i + 1)
            pylab.imshow(data[i,0,:,:], vmin = 0, vmax = 1)
            pylab.show()

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

    # apply the adjoint model (backward projection)
    imgs = am.backward(diff)

    print('---\n its backward projection norm: %e' % imgs.norm())

except error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
