import os
import pylab
import sys
import time

BUILD_PATH = os.environ.get('BUILD_PATH') + '/xGadgetron'
SRC_PATH = os.environ.get('SRC_PATH') + '/xGadgetron/pGadgetron'

sys.path.append(BUILD_PATH)
sys.path.append(SRC_PATH)

from pGadgetron import *
from pGadgets import *

try:
    # acquisitions will be read from this HDF file
    input_data = MR_Acquisitions('testdata.h5')
    #input_data = MR_Acquisitions('test_2D_2x.h5')

    # define gadgets
    gadget11 = Gadget('NoiseAdjustGadget')
    gadget12 = Gadget('AsymmetricEchoGadget')
    gadget13 = Gadget('RemoveROOversamplingGadget')
    gadget21 = Gadget('AcquisitionAccumulateTriggerGadget')
    gadget22 = Gadget('BucketToBufferGadget')
    gadget221 = Gadget('PrepRefGadget')
    #gadget23 = Gadget('SimpleReconGadget')
    gadget23 = Gadget('CartesianGrappaGadget')
    gadget241 = Gadget('FOVAdjustmentGadget')
    gadget242 = Gadget('ScalingGadget')
    gadget25 = Gadget('ImageArraySplitGadget')
    gadget31 = Gadget('ComplexToFloatGadget')
    gadget32 = Gadget('FloatToShortGadget')
    gadget6 = Gadget('ExtractGadget')

    acq_proc = AcquisitionsProcessor()
    acq_proc.add_gadget('g1', gadget11)
    acq_proc.add_gadget('g2', gadget12)
    acq_proc.add_gadget('g3', gadget13)
    print('processing acquisitions...')
    interim_data = acq_proc.process(input_data)

    # create reconstruction object
    recon = ImagesReconstructor()
    recon.add_gadget('g1', gadget21)
    recon.add_gadget('g2', gadget22)
    recon.add_gadget('g21', gadget221)
    recon.add_gadget('g3', gadget23)
    recon.add_gadget('g41', gadget241)
    recon.add_gadget('g42', gadget242)
    recon.add_gadget('g5', gadget25)
    #recon.add_gadget('g31', gadget31)
    #recon.add_gadget('g32', gadget32)
    recon.add_gadget('g6', gadget6)
    # connect to input data
    recon.set_input(interim_data)
    # perform reconstruction
    print('reconstructing...')
    recon.process()
    # get reconstructed images
    #interim_images = recon.get_output()
    images = recon.get_output()

    # build image post-processing chain
    img_proc0 = ImagesProcessor()
    images0 = img_proc0.process(images)
    #images0 = img_proc0.process(interim_images)

##    img_proc = ImagesProcessor()
##    img_proc.add_gadget('g6', gadget6)
##    #img_proc.add_gadget('g1', gadget31)
##    #img_proc.add_gadget('g2', gadget32)
##    # post-process reconstructed images
##    #interim_images.conversion_to_real(1)
##    print('processing images...')
##    images = img_proc.process(interim_images)

##    # plot obtained images
##    for i in range(images.number()):
##        data = images.image_as_array(i)
##        pylab.figure(i + 1)
##        pylab.imshow(data[0,0,:,:])
##        print('close the plot window to continue...')
##        pylab.show()

    nz = images.number()
    print('%d images' % nz)

    print('Please enter z-coordinate of the slice to view it')
    print('(a value outside the range [0 : %d] will stop this loop)'%(nz - 1))
    while True:
        s = str(input('z-coordinate: '))
        if len(s) < 1:
            break
        z = int(s)
        if z < 0 or z >= nz:
            break
        data = images.image_as_array(z)
        pylab.figure(z)
        pylab.imshow(data[0,0,:,:])
        print('delete the plot window to continue...')
        pylab.show()

except error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
