import math
import pylab
import sys
import time

sys.path.append('../../build/xGadgetron')
sys.path.append('../pGadgetron')
#import pGadgetron
#import pGadgets

from pGadgetron import *
from pGadgets import *

try:
    # acquisitions will be read from this HDF file
    input_data = MR_Acquisitions('testdata.h5')
    # use noiseless data to check the acquisition model
    # input_data = pGadgetron.ISMRMRDAcquisitions('ex_data.h5')

    # define gadgets
    gadget1 = RemoveROOversamplingGadget()
    gadget2 = SimpleReconstructionGadget()
    gadget3 = ExtractGadget()

    gadget2.set_property('trigger_dimension', 'repetition')
    gadget2.set_property('split_slices', 'true')

    acq_proc = AcquisitionsProcessor()
    acq_proc.add_gadget('g1', gadget1)
    print('processing acquisitions...')
    interim_data = acq_proc.process(input_data)

    # create reconstruction object
    recon = ImagesReconstructor()
    recon.add_gadget('g2', gadget2)
    # connect to input data
    recon.set_input(interim_data)
    # perform reconstruction
    print('reconstructing...')
    recon.process()
    # get reconstructed images
    interim_images = recon.get_output()

    # build image post-processing chain
    img_proc = ImagesProcessor()
    img_proc.add_gadget('g3', gadget3)
    # post-process reconstructed images
    print('processing images...')
    images = img_proc.process(interim_images)

    # plot obtained images
    for i in range(images.number()):
        data = interim_images.image_as_array(i)
        pylab.figure(i + 1)
        pylab.imshow(data[0,0,:,:])
        pylab.show()

except error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
