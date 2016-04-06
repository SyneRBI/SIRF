import pylab
import sys
import time

sys.path.append('../../build/xGadgetron')
sys.path.append('../pGadgetron')
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

    # build gadgets chain
    recon.add_gadget('g1', gadget1)
    recon.add_gadget('g2', gadget2)
    recon.add_gadget('g3', gadget3)
    recon.add_gadget('g4', gadget4)
    recon.add_gadget('g5', gadget5)
    recon.add_gadget('g6', gadget6)

    # connect to input data
    recon.set_input(input_data)
    # perform reconstruction
    recon.process()
    
    # get reconstructed images
    images = recon.get_output()

    # plot reconstructed images
    for i in range(images.number()):
        data = images.image_as_array(i)
        pylab.figure(i + 1)
        pylab.imshow(data[0,0,:,:])
        pylab.show()

    # write images to a new group in 'output1.h5'
    # named after the current date and time
    time_str = time.asctime()
    images.write('output1.h5', time_str)

except pGadgetron.error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
