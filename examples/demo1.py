'''
Hello World python client
'''

import sys
import time

sys.path.append('../../build/xGadgetron')
sys.path.append('../src')
import pGadgetron

time_str = time.asctime()                #date and time
print 'Hello World on', time_str

try:
	# acquisitions will be read from this HDF file
	input_data = pGadgetron.ISMRMRDataset('testdata.h5')
	
	# setup Gadgetron client connector
	conn = pGadgetron.ClientConnector()
	conn.set_timeout(10000)

	# reconstructions received from Gadgetron will be added to this HDF file 
	# as a new group named after the date and time of script execution
	conn.register_HDF_receiver('../../build/xGadgetron/output1.h5', time_str)
	
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

except pGadgetron.error as err:
    # display error information
    print 'Gadgetron exception occured:\n', err.value
