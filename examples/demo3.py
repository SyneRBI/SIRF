'''
Gadgetron python client
'''

import pylab
import sys
import time

sys.path.append('../../build/xGadgetron')
sys.path.append('../src')
import pGadgetron
import pGadgets

try:
	# acquisitions will be read from this HDF file
	input_data = pGadgetron.ISMRMRDataset('testdata.h5')
	
	# define gadgets
	reader = pGadgets.GadgetIsmrmrdAcquisitionMessageReader()
	writer = pGadgets.MRIImageWriter()
	gadget1 = pGadgets.RemoveROOversamplingGadget()
	gadget2 = pGadgets.AcquisitionAccumulateTriggerGadget()
	gadget3 = pGadgets.BucketToBufferGadget()
	gadget4 = pGadgets.SimpleReconGadget()
	gadget5 = pGadgets.ImageArraySplitGadget()
	gadget6 = pGadgets.ExtractGadget()
	gadget7 = pGadgets.ImageFinishGadget()

	# define gadgets chain
	gc = pGadgetron.GadgetChain()
	gc.addReader('r1', reader)
	gc.addWriter('w1', writer)
	gc.addGadget('g1', gadget1)
	gc.addGadget('g2', gadget2)
	gc.addGadget('g3', gadget3)
	gc.addGadget('g4', gadget4)
	gc.addGadget('g5', gadget5)
	gc.addGadget('g6', gadget6)
	gc.addGadget('g7', gadget7)

	# create rerconstruction object
	recon = pGadgetron.MRReconstructionDirect()
	# connect to input data
	recon.set_input(input_data)
	# perform reconstruction
	recon.process(gc)
	
	# get reconstructed images
	images = recon.get_output()

	# plot reconstructed images
	for i in range(images.number()):
		data = images.image_as_array(i)
		pylab.figure(i + 1)
		pylab.imshow(data[:,:,0])
		pylab.show()

	# write/append reconstructed images as a new data group in 'output3.h5' 
	# named after the current date and time
	#time_str = time.asctime()
	#images.write('../../build/xGadgetron/output3.h5', time_str)

except pGadgetron.error as err:
    # display error information
    print 'Gadgetron exception occured:\n', err.value
