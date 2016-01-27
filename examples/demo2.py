'''
Gadgetron python client
'''

import pylab
import sys
import time

sys.path.append('../../build/xGadgetron')
sys.path.append('../pGadgetron')
import pGadgetron

try:
	# acquisitions will be read from this HDF file
	input_data = pGadgetron.ISMRMRDataset('testdata.h5')
	
	# setup Gadgetron client connector
	conn = pGadgetron.ClientConnector()
	conn.set_timeout(10000)

	# reconstructions received from Gadgetron will be collected this list 
	img_list = pGadgetron.ImagesList()
	conn.register_images_receiver(img_list)
	
	# open Gadgetron connection
	conn.connect('localhost', '9002')

	# send gadget chain configuration defined by this xml file to Gadgetron
	conn.send_config_file('default.xml')
	
	# send acquisition parameters to Gadgetron
	conn.send_parameters(input_data.header)

	# start sending acquisition data to Gadgetron
	conn.send_acquisitions(input_data)
	
	# close connection to Gadgetron
	conn.disconnect()

	for i in range(img_list.size()):
		data = img_list.image_as_array(i)
		pylab.figure(i + 1)
		pylab.imshow(data[:,:,0])
		pylab.show()

	# write images to a new group in 'output2.h5' named after the current date and time
	time_str = time.asctime()
	img_list.write('../../build/xGadgetron/output2.h5', time_str)

except pGadgetron.error as err:
    # display error information
    print 'Gadgetron exception occured:\n', err.value
