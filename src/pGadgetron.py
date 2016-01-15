import numpy
import os
import time
import pygadgetron

class error(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def _check_status(handle):
    if pygadgetron.executionStatus(handle) != 0:
        msg = pygadgetron.executionError(handle)
        file = pygadgetron.executionErrorFile(handle)
        line = pygadgetron.executionErrorLine(handle)
        errorMsg = \
            repr(msg) + ' exception thrown at line ' + \
            repr(line) + ' of ' + file
        raise error(errorMsg)

class PyGadgetronObject:
	pass
#    def __init__(self):
#        self.handle = None
#    def __del__(self):
#        if self.handle is not None:
#            pygadgetron.deleteObject(self.handle)

class ImagesList(PyGadgetronObject):
	def __init__(self):
		self.handle = pygadgetron.newObject('ImagesList')
		_check_status(self.handle)
	def __del__(self):
		pygadgetron.deleteObject(self.handle)
	def size(self):
		return pygadgetron.cGT_numImages(self.handle)
	def number(self):
		return pygadgetron.cGT_numImages(self.handle)
	def write(self, out_file, out_group):
		conn = ClientConnector()
		handle = pygadgetron.cGT_writeImages\
			(self.handle, conn.handle, out_file, out_group)
		_check_status(handle)
	def image_as_array(self, im_num):
		dim = numpy.ndarray((3,), dtype = numpy.int32)
		pygadgetron.cGT_getImageDimensions(self.handle, im_num, dim.ctypes.data)
		nx = dim[0]
		ny = dim[1]
		nz = dim[2]
		if nx == 0 or ny == 0 or nz == 0:
			raise error('density data not available')
		array = numpy.ndarray((nx, ny, nz), dtype = numpy.float64)
		pygadgetron.cGT_getImageDataAsDoubleArray\
			(self.handle, im_num, array.ctypes.data)
		return array

class ClientConnector(PyGadgetronObject):
    def __init__(self):
        self.handle = None
        self.handle = pygadgetron.newObject('GTConnector')
        _check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)
    def set_timeout(self, timeout):
        handle = pygadgetron.cGT_setConnectionTimeout(self.handle, timeout)
        _check_status(handle)
        pygadgetron.deleteDataHandle(handle)
    def connect(self, host, port):
        handle = pygadgetron.cGT_connect(self.handle, host, port)
        _check_status(handle)
        pygadgetron.deleteDataHandle(handle)
    def disconnect(self):
        handle = pygadgetron.cGT_disconnect(self.handle)
        _check_status(handle)
        pygadgetron.deleteDataHandle(handle)
    def register_HDF_receiver(self, file, group):
        handle = pygadgetron.cGT_registerHDFReceiver(self.handle, file, group)
        _check_status(handle)
        pygadgetron.deleteDataHandle(handle)
    def register_images_receiver(self, imgs):
        handle = pygadgetron.cGT_registerImagesReceiver(self.handle, imgs.handle)
        _check_status(handle)
        pygadgetron.deleteDataHandle(handle)
    def send_config_file(self, file):
        handle = pygadgetron.cGT_sendConfigFile(self.handle, file)
        _check_status(handle)
        pygadgetron.deleteDataHandle(handle)
    def send_parameters(self, par):
        handle = pygadgetron.cGT_sendParameters(self.handle, par)
        _check_status(handle)
        pygadgetron.deleteDataHandle(handle)
    def send_acquisitions(self, acq):
        handle = pygadgetron.cGT_sendAcquisitions(self.handle, acq.handle)
        _check_status(handle)
        pygadgetron.deleteDataHandle(handle)
        
class ISMRMRDataset(PyGadgetronObject):
    def __init__(self, file):
        self.handle = None
        self.handle = pygadgetron.cGT_ISMRMRDatasetFromFile(file, '/dataset')
        _check_status(self.handle)
        self.header = pygadgetron.newObject('string')
        handle = pygadgetron.cGT_readISMRMRDatasetHeader(self.handle, self.header)
        _check_status(handle)
        pygadgetron.deleteDataHandle(handle)
    def __del__(self):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)
    def read_header(self):
        handle = pygadgetron.cGT_readISMRMRDatasetHeader(self.handle, self.header)
        _check_status(handle)
        pygadgetron.deleteDataHandle(handle)

class MRReconstruction(PyGadgetronObject):
	def __init__(self):
		self.conn = ClientConnector()
		self.images = ImagesList()
		self.conn.register_images_receiver(self.images)
		self.input_data = None
		self.gadgets_config = None
	def set_input(self, input_data):
		self.input_data = input_data
	def set_up(self, config_file):
		self.gadgets_config = config_file
	def process(self):
		if self.input_data is None:
			raise error('no input data')
		if self.gadgets_config is None:
			raise error('gadgets chain not defined')
		self.conn.connect('localhost', '9002')
		self.conn.send_config_file(self.gadgets_config)
		self.conn.send_parameters(self.input_data.header)
		self.conn.send_acquisitions(self.input_data)
		self.conn.disconnect()
	def get_output(self):
		return self.images

class MRReconstructionDirect(MRReconstruction):
	pass
