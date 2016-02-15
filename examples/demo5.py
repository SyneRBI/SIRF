import pylab
import sys
import time

sys.path.append('../../build/xGadgetron')
sys.path.append('../pGadgetron')
import pGadgetron
import pGadgets

try:
    # acquisitions will be read from this HDF file
    #print('creating acquisitions object...')
    input_data = pGadgetron.ISMRMRDAcquisitions('testdata.h5')

    # define gadgets
    gadget1 = pGadgets.RemoveROOversamplingGadget()
    gadget2 = pGadgets.SimpleReconstructionGadget()
    gadget3 = pGadgets.ExtractGadget()

    gadget2.set_property('trigger_dimension', 'repetition')
    gadget2.set_property('split_slices', 'true')

    #print('creating acquisitions processor object...')
    acq_proc = pGadgetron.AcquisitionsProcessor()
    print('scratch file:', acq_proc.acq_file)
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

    # plot obtained images
    for i in range(images.number()):
        data = images.image_as_array(i)
        pylab.figure(i + 1)
        pylab.imshow(data[:,:,0])
        pylab.show()

except pGadgetron.error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
