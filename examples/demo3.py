'''
Gadgetron python client
'''

import pylab
import sys
import time

sys.path.append('../../build/xGadgetron')
sys.path.append('../src')
import pGadgetron

try:
	# acquisitions will be read from this HDF file
	input_data = pGadgetron.ISMRMRDataset('testdata.h5')
	
	# create rerconstruction object
	recon = pGadgetron.MRReconstructionDirect()
	# connect to input data
	recon.set_input(input_data)
	# define gadgets chain
	recon.set_up('default.xml')
	# perform reconstruction
	recon.process()
	
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
	time_str = time.asctime()
	images.write('../../build/xGadgetron/output3.h5', time_str)

except pGadgetron.error as err:
    # display error information
    print 'Gadgetron exception occured:\n', err.value
