import pylab
import sys
import time

sys.path.append('../../build/xGadgetron')
sys.path.append('../pGadgetron')
import pGadgetron
import pGadgets

try:
    # acquisitions will be read from this HDF file
    input_data = pGadgetron.ISMRMRDAcquisitions('opismrmrd.h5')
    
    # define gadgets
    gadget1 = pGadgets.RemoveROOversamplingGadget()
    gadget2 = pGadgets.SimpleReconstructionGadget()
    gadget3 = pGadgets.ExtractGadget()

    # set gadgets parameters
    gadget2.set_property('trigger_dimension', 'repetition')
    gadget2.set_property('split_slices', 'true')

    # create reconstruction object
    recon = pGadgetron.ImagesReconstructor()

    # build gadgets chain
    recon.add_gadget('g1', gadget1)
    recon.add_gadget('g2', gadget2)
    recon.add_gadget('g3', gadget3)

    # connect to input data
    recon.set_input(input_data)
    # perform reconstruction
    recon.process()
    
    # get reconstructed images
    images = recon.get_output()

    # plot reconstructed images

    nz = images.number()
    print(nz, 'images')

    while True:
        s = str(input('enter z-coordinate: '))
        if len(s) < 1:
            break
        z = int(s)
        if z < 0 or z >= nz:
            break
        data = images.image_as_array(z)
        pylab.figure(z)
        pylab.imshow(data[:,:,0])
        pylab.show()

    # write images to a new group in 'output6.h5'
    # named after the current date and time
    time_str = time.asctime()
    images.write('output6.h5', time_str)

except pGadgetron.error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
