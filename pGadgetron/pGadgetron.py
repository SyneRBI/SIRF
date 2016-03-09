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
	
class ClientConnector(PyGadgetronObject):
    def __init__(self):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject('GTConnector')
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
        handle = pygadgetron.cGT_registerImagesReceiver\
            (self.handle, imgs.handle)
        _check_status(handle)
        pygadgetron.deleteDataHandle(handle)
    def config_gadget_chain(self, gc):
        handle = pygadgetron.cGT_configGadgetChain(self.handle, gc.handle)
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
    def send_images(self, img):
        handle = pygadgetron.cGT_sendImages(self.handle, img.handle)
        _check_status(handle)
        pygadgetron.deleteDataHandle(handle)
        
class ImagesContainer(PyGadgetronObject):
    def __init__(self):
        self.handle = None
##    def __init__(self, template = None):
##        self.handle = None
##        if template is None:
##            self.handle = pygadgetron.cGT_newObject('ImagesList')
##        else:
##            self.handle = pygadgetron.cGT_imagesCopy(template.handle)
##        _check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            #print('deleting images object...')
            pygadgetron.deleteObject(self.handle)
    def number(self):
        return pygadgetron.cGT_numImages(self.handle)
    def dot(self, images):
        handle = pygadgetron.cGT_imagesDot(self.handle, images.handle)
        _check_status(handle)
        re = pygadgetron.doubleReDataFromHandle(handle)
        im = pygadgetron.doubleImDataFromHandle(handle)
        return complex(re, im)
    def write(self, out_file, out_group):
        handle = pygadgetron.cGT_writeImages\
            (self.handle, out_file, out_group)
        _check_status(handle)
        pygadgetron.deleteDataHandle(handle)
    def image_as_array(self, im_num):
        dim = numpy.ndarray((4,), dtype = numpy.int32)
        pygadgetron.cGT_getImageDimensions\
            (self.handle, im_num, dim.ctypes.data)
        nx = dim[0]
        ny = dim[1]
        nz = dim[2]
        nc = dim[3]
        if nx == 0 or ny == 0 or nz == 0 or nc == 0:
            raise error('density data not available')
        array = numpy.ndarray((nx, ny, nz, nc), dtype = numpy.float64)
        pygadgetron.cGT_getImageDataAsDoubleArray\
            (self.handle, im_num, array.ctypes.data)
        return array

class AcquisitionsContainer(PyGadgetronObject):
    def __init__(self):
        self.handle = None
    def __del__(self):
        if self.handle is not None:
            #print('deleting acquisitions object...')
            pygadgetron.deleteObject(self.handle)
    def norm(self):
        handle = pygadgetron.cGT_acquisitionsNorm(self.handle)
        _check_status(handle)
        r = pygadgetron.doubleDataFromHandle(handle)
        pygadgetron.deleteDataHandle(handle)
        return r;
    def dot(self, acqs):
        handle = pygadgetron.cGT_acquisitionsDot(self.handle, acqs.handle)
        _check_status(handle)
        re = pygadgetron.doubleReDataFromHandle(handle)
        im = pygadgetron.doubleImDataFromHandle(handle)
        pygadgetron.deleteDataHandle(handle)
        return complex(re, im)

class ISMRMRDAcquisitions(AcquisitionsContainer):
    def __init__(self, file = None):
        self.handle = None
        if file is not None:
            self.handle = pygadgetron.cGT_ISMRMRDAcquisitionsFromFile(file)
            _check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)

def axpby(a, x, b, y):
    z = ISMRMRDAcquisitions()
    z.handle = pygadgetron.cGT_newAcquisitionsContainer(x.handle)
    handle = pygadgetron.cGT_acquisitionsAxpby\
        (a, x.handle, b, y.handle, z.handle)
    pygadgetron.deleteDataHandle(handle)
    return z;

class AcquisitionModel(PyGadgetronObject):
    def __init__(self, acqs, imgs):
        self.handle = None
        self.handle = \
            pygadgetron.cGT_AcquisitionModel(acqs.handle, imgs.handle)
        _check_status(self.handle)
        self.images = imgs
    def __del__(self):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)
    def forward(self, images):
        acqs = ISMRMRDAcquisitions()
        acqs.handle = pygadgetron.cGT_AcquisitionModelForward\
            (self.handle, images.handle)
        _check_status(acqs.handle)
        return acqs;
    def backward(self, acqs):
##        images = ImagesList(self.images)
##        handle = pygadgetron.cGT_AcquisitionModelBwd\
##            (self.handle, images.handle, acqs.handle)
##        _check_status(handle)
##        pygadgetron.deleteDataHandle(handle)
        images = ImagesContainer()
##        if images.handle is not None:
##            pygadgetron.deleteObject(images.handle)
        images.handle = pygadgetron.cGT_AcquisitionModelBackward\
            (self.handle, acqs.handle)
        _check_status(images.handle)
        return images

class GadgetChain(PyGadgetronObject):
    def __init__(self):
        self.handle = pygadgetron.cGT_newObject('GadgetChain')
        _check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)
    def add_reader(self, id, reader):
        handle = pygadgetron.cGT_addReader(self.handle, id, reader.handle)
        _check_status(handle)
        pygadgetron.deleteDataHandle(handle)
    def add_writer(self, id, writer):
        handle = pygadgetron.cGT_addWriter(self.handle, id, writer.handle)
        _check_status(handle)
        pygadgetron.deleteDataHandle(handle)
    def add_gadget(self, id, gadget):
        handle = pygadgetron.cGT_addGadget(self.handle, id, gadget.handle)
        _check_status(handle)
        pygadgetron.deleteDataHandle(handle)

class ImagesReconstructor(GadgetChain):
    def __init__(self):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject('ImageReconstructor')
        _check_status(self.handle)
        self.input_data = None
    def __del__(self):
        if self.handle is not None:
            #print('deleting reconstructor object...')
            pygadgetron.deleteObject(self.handle)
    def set_input(self, input_data):
        self.input_data = input_data
    def process(self):
        if self.input_data is None:
            raise error('no input data')
        handle = pygadgetron.cGT_reconstructImages\
             (self.handle, self.input_data.handle)
        _check_status(handle)
        pygadgetron.deleteDataHandle(handle)
    def get_output(self):
        images = ImagesContainer()
##        if images.handle is not None:
##            pygadgetron.deleteObject(images.handle)
        images.handle = pygadgetron.cGT_reconstructedImages(self.handle)
        _check_status(images.handle)
        return images

class ImagesProcessor(GadgetChain):
    def __init__(self):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject('ImagesProcessor')
        _check_status(self.handle)
        self.input_data = None
    def __del__(self):
        if self.handle is not None:
            #print('deleting image processor object...')
            pygadgetron.deleteObject(self.handle)
##    def set_input(self, input_data):
##        self.input_data = input_data
    def process(self, input_data):
##        if self.input_data is None:
##            raise error('no input data')
        images = ImagesContainer()
##        if images.handle is not None:
##            pygadgetron.deleteObject(images.handle)
        images.handle = pygadgetron.cGT_processImages\
             (self.handle, input_data.handle)
        _check_status(images.handle)
        return images
##    def get_output(self):
##        images = ImagesList()
##        if images.handle is not None:
##            pygadgetron.deleteObject(images.handle)
##        images.handle = pygadgetron.cGT_reconstructedImagesList(self.handle)
##        _check_status(images.handle)
##        return images

class AcquisitionsProcessor(GadgetChain):
    def __init__(self):
        self.handle = None
        self.handle = pygadgetron.cGT_acquisitionsProcessor()
        _check_status(self.handle)
        self.input_data = None
    def __del__(self):
        if self.handle is not None:
            #print('deleting acquisitions processor object...')
            pygadgetron.deleteObject(self.handle)
    def process(self, input_data):
        acquisitions = ISMRMRDAcquisitions()
        acquisitions.handle = pygadgetron.cGT_processAcquisitions\
             (self.handle, input_data.handle)
        _check_status(acquisitions.handle)
        return acquisitions

    
