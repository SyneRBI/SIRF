import os
import pylab
import sys
import time

BUILD_PATH = os.environ.get('BUILD_PATH') + '/xGadgetron'
SRC_PATH = os.environ.get('SRC_PATH') + '/xGadgetron/pGadgetron'

sys.path.append(BUILD_PATH)
sys.path.append(SRC_PATH)

import pGadgetron
import pGadgets

try:
    # acquisitions will be read from this HDF file
    input_data = pGadgetron.MR_Acquisitions('testdata.h5')
    
    # define gadgets
    gadget1 = pGadgets.RemoveROOversamplingGadget()
    gadget2 = pGadgets.AcquisitionAccumulateTriggerGadget()
    gadget3 = pGadgets.BucketToBufferGadget()
    gadget4 = pGadgets.SimpleReconGadget()
    gadget5 = pGadgets.ImageArraySplitGadget()
    gadget6 = pGadgets.ExtractGadget()

    gadget2.set_property('trigger_dimension', 'repetition')
    gadget3.set_property('split_slices', 'true')

    # create reconstruction object
    recon = pGadgetron.ImagesReconstructor()

    # build reconstruction chain
    recon.add_gadget('g1', gadget1)
    recon.add_gadget('g2', gadget2)
    recon.add_gadget('g3', gadget3)
    recon.add_gadget('g4', gadget4)
    recon.add_gadget('g5', gadget5)

    # connect to input data
    recon.set_input(input_data)
    # perform reconstruction
    recon.process()

    # get reconstructed images
    imgs = recon.get_output()

    # build image processing chain
    proc = pGadgetron.ImagesProcessor()
    proc.add_gadget('g6', gadget6)

    images = proc.process(imgs)

    # plot reconstructed images
    for i in range(images.number()):
        data = images.image_as_array(i)
        pylab.figure(i + 1)
        pylab.imshow(data[0,0,:,:])
        print('delete the plot window to continue...')
        pylab.show()

except pGadgetron.error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
