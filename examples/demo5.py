import math
import pylab
import sys
import time

sys.path.append('../../build/xGadgetron')
sys.path.append('../pGadgetron')
import pGadgetron
import pGadgets

try:
    # acquisitions will be read from this HDF file
    input_data = pGadgetron.ISMRMRDAcquisitions('testdata.h5')
    # use noiseless data to check the acquisition model
    # input_data = pGadgetron.ISMRMRDAcquisitions('ex_data.h5')

    # define gadgets
    gadget1 = pGadgets.RemoveROOversamplingGadget()
    gadget2 = pGadgets.SimpleReconstructionGadget()
    gadget3 = pGadgets.ExtractGadget()

    gadget2.set_property('trigger_dimension', 'repetition')
    gadget2.set_property('split_slices', 'true')

    acq_proc = pGadgetron.AcquisitionsProcessor()
    acq_proc.add_gadget('g1', gadget1)
    print('processing acquisitions...')
    interim_data = acq_proc.process(input_data)

    # create reconstruction object
    recon = pGadgetron.ImagesReconstructor()
    recon.add_gadget('g2', gadget2)
    # connect to input data
    recon.set_input(interim_data)
    # perform reconstruction
    print('reconstructing...')
    recon.process()
    # get reconstructed images
    interim_images = recon.get_output()

    # build image post-processing chain
    img_proc = pGadgetron.ImagesProcessor()
    img_proc.add_gadget('g3', gadget3)
    # post-process reconstructed images
    print('processing images...')
    images = img_proc.process(interim_images)

    # create acquisition model based on the acquisition parameters
    # stored in input_data and image parameters stored in interim_images
    am = pGadgetron.AcquisitionModel(input_data, interim_images)

    # use the acquisition model (forward projection) to produce acquisitions
    acqs = am.forward(interim_images)

    # compute the difference between real and modelled acquisitions
    a = -acqs.dot(input_data) / input_data.dot(input_data)
    a = a.real
    b = 1.0
    diff = pGadgetron.axpby(a, input_data, b, acqs)
    print(diff.norm()/acqs.norm())

    # apply the adjoint model (backward projection)
    imgs = am.backward(diff)

    # test that the backward projection is the adjoint of forward
    # on x = diff and y = interim_images
    print('(x, F y) = ', diff.dot(acqs))
    print('(B x, y) = ', imgs.dot(interim_images))

    # plot obtained images
    for i in range(images.number()):
        data = images.image_as_array(i)
        pylab.figure(i + 1)
        pylab.imshow(data[:,:,0,0])
        pylab.show()

except pGadgetron.error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
