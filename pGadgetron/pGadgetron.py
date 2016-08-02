import numpy
import os
try:
    import pylab
    HAVE_PYLAB = True
except:
    HAVE_PYLAB = False
import time

import pygadgetron

ISMRMRD_IMTYPE_MAGNITUDE = 1
ISMRMRD_IMTYPE_PHASE     = 2
ISMRMRD_IMTYPE_REAL      = 3
ISMRMRD_IMTYPE_IMAG      = 4

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

def _setParameter(hs, set, par, hv):
    h = pygadgetron.cGT_setParameter(hs, set, par, hv)
    _check_status(h)
    pygadgetron.deleteDataHandle(h)
def _set_int_par(handle, set, par, value):
    h = pygadgetron.intDataHandle(value)
    _setParameter(handle, set, par, h)
    pygadgetron.deleteDataHandle(h)
def _int_par(handle, set, par):
    h = pygadgetron.cGT_parameter(handle, set, par)
    _check_status(h)
    value = pygadgetron.intDataFromHandle(h)
    pygadgetron.deleteDataHandle(h)
    return value

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
        
class DataContainer(PyGadgetronObject):
    def __init__(self):
        self.handle = None
        self.is_real = False
    def __del__(self):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)
    def number(self):
        handle = pygadgetron.cGT_dataItems(self.handle)
        _check_status(handle)
        n = pygadgetron.intDataFromHandle(handle)
        pygadgetron.deleteDataHandle(handle)
        return n
    def real(self):
        if self.is_real:
            return self
        handle = pygadgetron.cGT_newObject('ExtractRealImagesProcessor')
        _check_status(handle)
        real_images = ImagesContainer()
        real_images.handle = pygadgetron.cGT_processImages\
             (handle, self.handle)
        _check_status(real_images.handle)
        pygadgetron.deleteObject(handle)
        real_images.is_real = True
        return real_images
    def norm(self):
        handle = pygadgetron.cGT_norm(self.handle)
        _check_status(handle)
        r = pygadgetron.doubleDataFromHandle(handle)
        pygadgetron.deleteDataHandle(handle)
        return r;
    def dot(self, acqs):
        handle = pygadgetron.cGT_dot(self.handle, acqs.handle)
        _check_status(handle)
        re = pygadgetron.doubleReDataFromHandle(handle)
        im = pygadgetron.doubleImDataFromHandle(handle)
        pygadgetron.deleteDataHandle(handle)
        return complex(re, im)
    def __add__(self, other):
        z = DataContainer()
        z.handle = pygadgetron.cGT_axpby\
            (1.0, 0.0, self.handle, 1.0, 0.0, other.handle)
        return z;
    def __sub__(self, other):
        z = DataContainer()
        z.handle = pygadgetron.cGT_axpby\
            (1.0, 0.0, self.handle, -1.0, 0.0, other.handle)
        return z;
    def __mul__(self, other):
        #print(type(other))
        if isinstance(other, DataContainer):
            return self.dot(other)
        elif type(other) == type(complex(0,0)):
            z = DataContainer()
            z.handle = pygadgetron.cGT_axpby\
                (other.real, other.imag, self.handle, 0, 0, self.handle)
            return z;
        elif type(other) == type(0.0):
            z = DataContainer()
            z.handle = pygadgetron.cGT_axpby\
                (other, 0, self.handle, 0, 0, self.handle)
            return z;
        else:
            raise error('wrong multiplier')
    def __rmul__(self, other):
        if type(other) == type(complex(0,0)):
            z = DataContainer()
            z.handle = pygadgetron.cGT_axpby\
                (other.real, other.imag, self.handle, 0, 0, self.handle)
            return z;
        elif type(other) == type(0.0):
            z = DataContainer()
            z.handle = pygadgetron.cGT_axpby\
                (other, 0, self.handle, 0, 0, self.handle)
            return z;
        else:
            raise error('wrong multiplier')
    @staticmethod
    def axpby(a, x, b, y):
        z = DataContainer()
        z.handle = pygadgetron.cGT_axpby\
            (a.real, a.imag, x.handle, b.real, b.imag, y.handle)
        return z;

class MR_CoilSensitivityMaps(DataContainer):
    def __init__(self):
        self.handle = None
        self.smoothness = 0
    def __del__(self):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)
    def read(self, file):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)
        self.handle = pygadgetron.cGT_CoilSensitivities(file)
        _check_status(self.handle)
    def set_smoothness(self, s):
        self.smoothness = s
    def calculate(self, acqs):
        if acqs.is_sorted() is False:
            print('WARNING: acquisitions may be in a wrong order')
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)
        self.handle = pygadgetron.cGT_CoilSensitivities('')
        _check_status(self.handle)
        _set_int_par\
            (self.handle, 'coil_sensitivity', 'smoothness', self.smoothness)
        handle = pygadgetron.cGT_computeCoilSensitivities\
            (self.handle, acqs.handle)
        _check_status(handle)
        pygadgetron.deleteDataHandle(handle)
    def csm_as_array(self, csm_num):
        dim = numpy.ndarray((4,), dtype = numpy.int32)
        pygadgetron.cGT_getCSMDimensions\
            (self.handle, csm_num, dim.ctypes.data)
        nx = dim[0]
        ny = dim[1]
        nz = dim[2]
        nc = dim[3]
        if nx == 0 or ny == 0 or nz == 0 or nc == 0:
            raise error('image data not available')
        array = numpy.ndarray((nc, nz, ny, nx), dtype = numpy.float64)
        pygadgetron.cGT_getCSMDataAbs\
            (self.handle, csm_num, array.ctypes.data)
        return array
    def csm_as_arrays(self, csm_num):
        dim = numpy.ndarray((4,), dtype = numpy.int32)
        pygadgetron.cGT_getCSMDimensions\
            (self.handle, csm_num, dim.ctypes.data)
        nx = dim[0]
        ny = dim[1]
        nz = dim[2]
        nc = dim[3]
        if nx == 0 or ny == 0 or nz == 0 or nc == 0:
            raise error('image data not available')
        re = numpy.ndarray((nc, nz, ny, nx), dtype = numpy.float64)
        im = numpy.ndarray((nc, nz, ny, nx), dtype = numpy.float64)
        pygadgetron.cGT_getCSMData\
            (self.handle, csm_num, re.ctypes.data, im.ctypes.data)
        return re, im

class ImagesContainer(DataContainer):
    def __init__(self):
        self.handle = None
        self.is_real = False
    def __del__(self):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)
    def types(self):
        handle = pygadgetron.cGT_imageTypes(self.handle)
        _check_status(handle)
        n = pygadgetron.intDataFromHandle(handle)
        pygadgetron.deleteDataHandle(handle)
        return n
    def show(self):
        if not HAVE_PYLAB:
            print('pylab not found')
            return
        ni = self.number()
        print('%d images' % ni)
        print('Please enter the number of the image to view')
        print('(a value outside the range [1 : %d] will stop this loop)' % ni)
        while True:
            s = str(input('image: '))
            if len(s) < 1:
                break
            i = int(s)
            if i < 1 or i > ni:
                break
            data = self.image_as_array(i - 1)
            pylab.figure(i)
            pylab.title('image %d' % i)
            pylab.imshow(data[0,0,:,:])
            print('Close Figure %d window to continue...' % i)
            pylab.show()
    def write(self, out_file, out_group):
        handle = pygadgetron.cGT_writeImages\
            (self.handle, out_file, out_group)
        _check_status(handle)
        pygadgetron.deleteDataHandle(handle)
    def select(self, inc = 1, off = 0):
        images = ImagesContainer()
        images.handle = pygadgetron.cGT_selectImages(self.handle, inc, off)
        _check_status(images.handle)
        return images
    def conversion_to_real(self, ctype):
        pygadgetron.cGT_setImageToRealConversion(self.handle, ctype)
    def image_as_array(self, im_num):
        dim = numpy.ndarray((4,), dtype = numpy.int32)
        pygadgetron.cGT_getImageDimensions\
            (self.handle, im_num, dim.ctypes.data)
        nx = dim[0]
        ny = dim[1]
        nz = dim[2]
        nc = dim[3]
        if nx == 0 or ny == 0 or nz == 0 or nc == 0:
            raise error('image data not available')
##        array = numpy.ndarray((nx, ny, nz, nc), dtype = numpy.float64)
        array = numpy.ndarray((nc, nz, ny, nx), dtype = numpy.float64)
        pygadgetron.cGT_getImageDataAsDoubleArray\
            (self.handle, im_num, array.ctypes.data)
        return array
    def image_as_complex_array(self, im_num):
        dim = numpy.ndarray((4,), dtype = numpy.int32)
        pygadgetron.cGT_getImageDimensions\
            (self.handle, im_num, dim.ctypes.data)
        nx = dim[0]
        ny = dim[1]
        nz = dim[2]
        nc = dim[3]
        if nx == 0 or ny == 0 or nz == 0 or nc == 0:
            raise error('image data not available')
        array = numpy.ndarray((nc, nz, ny, nx), dtype = numpy.complex64)
        pygadgetron.cGT_getImageDataAsComplexArray\
            (self.handle, im_num, array.ctypes.data)
        return array

class AcquisitionsContainer(DataContainer):
    def __init__(self):
        self.handle = None
        self.sorted = False
    def __del__(self):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)
    def sort(self):
        handle = pygadgetron.cGT_orderAcquisitions(self.handle)
        _check_status(handle)
        pygadgetron.deleteDataHandle(handle)
        self.sorted = True
    def is_sorted(self):
        return self.sorted
    def process(self, list):
        ap = AcquisitionsProcessor(list)
        return ap.process(self)
    def acquisition(self, num):
        acq = MR_Acquisition()
        acq.handle = pygadgetron.cGT_acquisitionFromContainer(self.handle, num)
        return acq
    def slice_dimensions(self):
        dim = numpy.ndarray((3,), dtype = numpy.int32)
        pygadgetron.cGT_getAcquisitionsDimensions(self.handle, dim.ctypes.data)
        ns = dim[0]
        ny = dim[1]
        nc = dim[2]
        return ns, ny, nc
    def slice_as_array(self, num):
        dim = numpy.ndarray((3,), dtype = numpy.int32)
        pygadgetron.cGT_getAcquisitionsDimensions(self.handle, dim.ctypes.data)
        ns = dim[0]
        ny = dim[1]
        nc = dim[2]
        #print( ns, ny, nc )
        re = numpy.ndarray((nc, ny, ns), dtype = numpy.float64)
        im = numpy.ndarray((nc, ny, ns), dtype = numpy.float64)
        pygadgetron.cGT_getAcquisitionsData\
            (self.handle, num, re.ctypes.data, im.ctypes.data)
        return re + 1j*im

class MR_Acquisition(PyGadgetronObject):
    def __init__(self, file = None):
        self.handle = None
    def __del__(self):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)
    def flags(self):
        return _int_par(self.handle, 'acquisition', 'flags')
    def number_of_samples(self):
        return _int_par(self.handle, 'acquisition', 'number_of_samples')
    def active_channels(self):
        return _int_par(self.handle, 'acquisition', 'active_channels')
    def trajectory_dimensions(self):
        return _int_par(self.handle, 'acquisition', 'trajectory_dimensions')
    def idx_kspace_encode_step_1(self):
        return _int_par(self.handle, 'acquisition', 'idx_kspace_encode_step_1')
    def idx_repetition(self):
        return _int_par(self.handle, 'acquisition', 'idx_repetition')
    def idx_slice(self):
        return _int_par(self.handle, 'acquisition', 'idx_slice')

class MR_Acquisitions(AcquisitionsContainer):
    def __init__(self, file = None):
        self.handle = None
        self.sorted = False
        if file is not None:
            self.handle = pygadgetron.cGT_ISMRMRDAcquisitionsFromFile(file)
            _check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)

class MR_AcquisitionModel(PyGadgetronObject):
    def __init__(self, acqs, imgs):
        self.handle = None
        self.handle = \
            pygadgetron.cGT_AcquisitionModel(acqs.handle, imgs.handle)
        _check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)
    def set_coil_sensitivity_maps(self, csm):
        handle = pygadgetron.cGT_setCSMs(self.handle, csm.handle)
        _check_status(handle)
        pygadgetron.deleteDataHandle(handle)
        self.csms = 'set'
    def forward(self, images):
        acqs = MR_Acquisitions()
        acqs.handle = pygadgetron.cGT_AcquisitionModelForward\
            (self.handle, images.handle)
        _check_status(acqs.handle)
        return acqs;
    def backward(self, acqs):
        images = ImagesContainer()
        images.handle = pygadgetron.cGT_AcquisitionModelBackward\
            (self.handle, acqs.handle)
        _check_status(images.handle)
        return images

class Gadget(PyGadgetronObject):
    def __init__(self, name):
        self.handle = None
        i = name.find('(')
        if i > -1:
            j = name.find(')', i)
            prop = name[i + 1 : j]
            name = name[: i]
            #name.rstrip() # does not work, have to do this way:
            i = name.find(' ')
            if i > -1:
                name = name[:i]
            i = 0
            #print(name, prop)
        self.handle = pygadgetron.cGT_newObject(name)
        _check_status(self.handle)
        if i > -1:
            self.set_properties(prop)
    def __del__(self):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)
    def set_property(self, prop, value):
        handle = pygadgetron.cGT_setGadgetProperty(self.handle, prop, value)
        _check_status(handle)
        pygadgetron.deleteDataHandle(handle)
    def set_properties(self, prop):
        handle = pygadgetron.cGT_setGadgetProperties(self.handle, prop)
        _check_status(handle)
        pygadgetron.deleteDataHandle(handle)

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
    def __init__(self, list = None):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject('ImagesReconstructor')
        _check_status(self.handle)
        self.input_data = None
        if list is None:
            return
        for i in range(len(list)):
            self.add_gadget('g' + repr(i + 1), Gadget(list[i]))
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
        images.handle = pygadgetron.cGT_reconstructedImages(self.handle)
        _check_status(images.handle)
        return images

class ImagesProcessor(GadgetChain):
    def __init__(self, list = None):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject('ImagesProcessor')
        _check_status(self.handle)
        self.input_data = None
        if list is None:
            return
        for i in range(len(list)):
            self.add_gadget('g' + repr(i + 1), Gadget(list[i]))
    def __del__(self):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)
##    def set_input(self, input_data):
##        self.input_data = input_data
    def process(self, input_data):
##        if self.input_data is None:
##            raise error('no input data')
        images = ImagesContainer()
        images.handle = pygadgetron.cGT_processImages\
             (self.handle, input_data.handle)
        _check_status(images.handle)
        return images
##    def get_output(self):
##        images = ImagesContainer()
##        images.handle = pygadgetron.cGT_reconstructedImagesList(self.handle)
##        _check_status(images.handle)
##        return images

class AcquisitionsProcessor(GadgetChain):
    def __init__(self, list = None):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject('AcquisitionsProcessor')
        _check_status(self.handle)
        self.input_data = None
        if list is None:
            return
        for i in range(len(list)):
            self.add_gadget('g' + repr(i + 1), Gadget(list[i]))
    def __del__(self):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)
    def process(self, input_data):
        acquisitions = AcquisitionsContainer()
        acquisitions.handle = pygadgetron.cGT_processAcquisitions\
             (self.handle, input_data.handle)
        _check_status(acquisitions.handle)
        return acquisitions

class MR_BasicReconstruction(ImagesReconstructor):
    def __init__(self):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject('SimpleReconstructionProcessor')
        _check_status(self.handle)
        self.input_data = None
    def __del__(self):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)
    
class MR_BasicGRAPPAReconstruction(ImagesReconstructor):
    def __init__(self):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject\
            ('SimpleGRAPPAReconstructionProcessor')
        _check_status(self.handle)
        self.input_data = None
    def __del__(self):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)
    def get_output(self):
        output = ImagesReconstructor.get_output(self)
        images = output.select(2)
        gfactors = output.select(2,1)
        return images, gfactors
    
def MR_remove_x_oversampling(input_data):
    handle = pygadgetron.cGT_newObject('RemoveOversamplingProcessor')
    _check_status(handle)
    output_data = AcquisitionsContainer()
    output_data.handle = pygadgetron.cGT_processAcquisitions\
         (handle, input_data.handle)
    _check_status(output_data.handle)
    pygadgetron.deleteObject(handle)
    return output_data

def MR_extract_real_images(complex_images):
    handle = pygadgetron.cGT_newObject('ExtractRealImagesProcessor')
    _check_status(handle)
    real_images = ImagesContainer()
    real_images.handle = pygadgetron.cGT_processImages\
         (handle, complex_images.handle)
    _check_status(real_images.handle)
    pygadgetron.deleteObject(handle)
    return real_images

