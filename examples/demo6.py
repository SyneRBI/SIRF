'''
Lower level single-chain GRAPPA reconstruction demo.
'''

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

    # define gadgets
    gadget11 = Gadget('NoiseAdjustGadget')
    gadget12 = Gadget('AsymmetricEchoGadget')
    gadget13 = Gadget('RemoveROOversamplingGadget')
    gadget21 = Gadget('AcquisitionAccumulateTriggerGadget')
    gadget22 = Gadget('BucketToBufferGadget')
    gadget23 = Gadget('PrepRefGadget')
    gadget24 = Gadget('CartesianGrappaGadget')
    gadget25 = Gadget('FOVAdjustmentGadget')
    gadget26 = Gadget('ScalingGadget')
    gadget27 = Gadget('ImageArraySplitGadget')
    gadget31 = Gadget('ComplexToFloatGadget')
    gadget32 = Gadget('FloatToShortGadget')

    acq_proc = AcquisitionsProcessor()
    acq_proc.add_gadget('g1', gadget11)
    acq_proc.add_gadget('g2', gadget12)
    acq_proc.add_gadget('g3', gadget13)
    print('pre-processing acquisitions...')
    preprocessed_data = acq_proc.process(input_data)

    # create reconstruction object
    recon = ImagesReconstructor()
    recon.add_gadget('g1', gadget21)
    recon.add_gadget('g2', gadget22)
    recon.add_gadget('g3', gadget23)
    recon.add_gadget('g4', gadget24)
    recon.add_gadget('g5', gadget25)
    recon.add_gadget('g6', gadget26)
    recon.add_gadget('g7', gadget27)
    # connect to input data
    recon.set_input(preprocessed_data)
    # perform reconstruction
    print('reconstructing...')
    recon.process()
    # get reconstructed images
    complex_images = recon.get_output()

    img_proc = ImagesProcessor()
    img_proc.add_gadget('g1', gadget31)
    img_proc.add_gadget('g2', gadget32)
    # post-process reconstructed images
    complex_images.conversion_to_real(1)
    print('processing images...')
    images = img_proc.process(complex_images)

    nz = images.number()
    print('%d images reconstructed.' % nz)

    print('Enter z-coordinate of the slice to view it')
    print('(a value outside the range [0 : %d] will stop this loop)'%(nz - 1))
    while True:
        s = str(input('z-coordinate: '))
        if len(s) < 1:
            break
        z = int(s)
        if z < 0 or z >= nz:
            break
        data = images.image_as_array(z)
        pylab.figure(z + 1)
        pylab.imshow(data[0,0,:,:])
        print('Close Figure %d window to continue...' % (z + 1))
        pylab.show()

    # write images to a new group in 'output10.h5'
    # named after the current date and time
    print('appending output10.h5...')
    time_str = time.asctime()
    images.write('output10.h5', time_str)

except error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
