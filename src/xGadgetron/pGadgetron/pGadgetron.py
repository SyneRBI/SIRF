''' 
Object-Oriented wrap for the cGadgetron-to-Python interface pygadgetron.py
'''
import numpy
import os
try:
    import pylab
    HAVE_PYLAB = True
except:
    HAVE_PYLAB = False
import sys
import time
try:
    from ismrmrdtools import coils
    HAVE_ISMRMRDTOOLS = True
except:
    print('ismrmrd-python-tools not installed')
    HAVE_ISMRMRDTOOLS = False

from pUtil import *
import pyiutil
import pygadgetron

# max number of acquisitions dimensiona
MAX_ACQ_DIMENSIONS = 16

# mask for image-related acquisitions
IMAGE_DATA_MASK = 0x11BFFFF

# image type
ISMRMRD_IMTYPE_MAGNITUDE = 1
ISMRMRD_IMTYPE_PHASE     = 2
ISMRMRD_IMTYPE_REAL      = 3
ISMRMRD_IMTYPE_IMAG      = 4

#image data type
ISMRMRD_USHORT   = 1 ##, /**< corresponds to uint16_t */
ISMRMRD_SHORT    = 2 ##, /**< corresponds to int16_t */
ISMRMRD_UINT     = 3 ##, /**< corresponds to uint32_t */
ISMRMRD_INT      = 4 ##, /**< corresponds to int32_t */
ISMRMRD_FLOAT    = 5 ##, /**< corresponds to float */
ISMRMRD_DOUBLE   = 6 ##, /**< corresponds to double */
ISMRMRD_CXFLOAT  = 7 ##, /**< corresponds to complex float */
ISMRMRD_CXDOUBLE = 8 ##  /**< corresponds to complex double */

###########################################################
############ Utilities for internal use only ##############
def _setParameter(hs, set, par, hv):
    h = pygadgetron.cGT_setParameter(hs, set, par, hv)
    check_status(h)
    pyiutil.deleteDataHandle(h)
def _set_int_par(handle, set, par, value):
    h = pyiutil.intDataHandle(value)
    _setParameter(handle, set, par, h)
    pyiutil.deleteDataHandle(h)
def _int_par(handle, set, par):
    h = pygadgetron.cGT_parameter(handle, set, par)
    check_status(h)
    value = pyiutil.intDataFromHandle(h)
    pyiutil.deleteDataHandle(h)
    return value
def _char_par(handle, set, par):
    h = pygadgetron.cGT_parameter(handle, set, par)
    check_status(h)
    value = pyiutil.charDataFromHandle(h)
    pyiutil.deleteDataHandle(h)
    return value
def _parameterHandle(hs, set, par):
    handle = pygadgetron.cGT_parameter(hs, set, par)
    check_status(handle)
    return handle
###########################################################

# data path finding helper functions
def mr_data_path():
    '''
    Returns default path to MR raw data files.
    '''
    return petmr_data_path('mr')
def raw_data_path():
    '''
    Returns default path to MR raw data files.
    '''
    return petmr_data_path('mr')

# abstract base class for everything, just in case we need some functionality
# common to everything
class PyGadgetronObject:
    pass

# low-level client functionality
# likely to be obsolete- not used for a long time
class ClientConnector(PyGadgetronObject):
    def __init__(self):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject('GTConnector')
        check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteObject(self.handle)
    def set_timeout(self, timeout):
        handle = pygadgetron.cGT_setConnectionTimeout(self.handle, timeout)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
    def connect(self, host, port):
        handle = pygadgetron.cGT_connect(self.handle, host, port)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
    def disconnect(self):
        handle = pygadgetron.cGT_disconnect(self.handle)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
    def register_images_receiver(self, imgs):
        handle = pygadgetron.cGT_registerImagesReceiver\
            (self.handle, imgs.handle)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
    def config_gadget_chain(self, gc):
        handle = pygadgetron.cGT_configGadgetChain(self.handle, gc.handle)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
    def send_config_file(self, file):
        handle = pygadgetron.cGT_sendConfigFile(self.handle, file)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
    def send_parameters(self, par):
        handle = pygadgetron.cGT_sendParameters(self.handle, par)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
    def send_acquisitions(self, acq):
        handle = pygadgetron.cGT_sendAcquisitions(self.handle, acq.handle)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
    def send_images(self, img):
        handle = pygadgetron.cGT_sendImages(self.handle, img.handle)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)

# base class for data container classes
# must actually be moved together with the respective C interface
# from xGadgetron/cGadgetron to common
class DataContainer(PyGadgetronObject):
    '''
    Class for an abstract data container.
    '''
    def __init__(self):
        self.handle = None
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteObject(self.handle)
    def same_object(self):
        '''
        Returns an object of the same type as self.
        '''
        # a kludge replacing purely abstract virtual method
        # necessitated by different syntax for such methods in Python 2 and 3
        if isinstance(self, ImageData):
            return ImageData()
        if isinstance(self, AcquisitionData):
            return AcquisitionData()
        return DataContainer()
    def number(self):
        '''
        Returns the number of items in the container.
        '''
        handle = pygadgetron.cGT_dataItems(self.handle)
        check_status(handle)
        n = pyiutil.intDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return n
    def norm(self):
        '''
        Returns the 2-norm of the container data viewed as a vector.
        '''
        handle = pygadgetron.cGT_norm(self.handle)
        check_status(handle)
        r = pyiutil.doubleDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return r;
    def dot(self, other):
        '''
        Returns the dot product of the container data with another container 
        data viewed as vectors.
        other: DataContainer
        '''
        handle = pygadgetron.cGT_dot(self.handle, other.handle)
        check_status(handle)
        re = pyiutil.doubleReDataFromHandle(handle)
        im = pyiutil.doubleImDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return complex(re, im)
    def __add__(self, other):
        '''
        Overloads + for data containers.
        Returns the sum of the container data with another container 
        data viewed as vectors.
        other: DataContainer
        '''
        z = self.same_object()
        z.handle = pygadgetron.cGT_axpby\
            (1.0, 0.0, self.handle, 1.0, 0.0, other.handle)
        return z;
    def __sub__(self, other):
        '''
        Overloads - for data containers.
        Returns the difference of the container data with another container 
        data viewed as vectors.
        other: DataContainer
        '''
        z = self.same_object()
        z.handle = pygadgetron.cGT_axpby\
            (1.0, 0.0, self.handle, -1.0, 0.0, other.handle)
        return z;
    def __mul__(self, other):
        '''
        Overloads * for data containers multiplication by a scalar.
        Returns the product self*s where s is a scalar factor.
        other: a real or complex scalar
        '''
        if isinstance(other, DataContainer):
            return self.dot(other)
        z = self.same_object()
        if type(other) == type(complex(0,0)):
            z.handle = pygadgetron.cGT_axpby\
                (other.real, other.imag, self.handle, 0, 0, self.handle)
            return z;
        elif type(other) == type(0.0):
            z.handle = pygadgetron.cGT_axpby\
                (other, 0, self.handle, 0, 0, self.handle)
            return z;
        else:
            raise error('wrong multiplier')
    def __rmul__(self, other):
        '''
        Overloads * for data containers multiplication by a scalar.
        Returns the product s*self where s is a scalar factor.
        other: a real or complex scalar
        '''
        z = self.same_object()
        if type(other) == type(complex(0,0)):
            z.handle = pygadgetron.cGT_axpby\
                (other.real, other.imag, self.handle, 0, 0, self.handle)
            return z;
        elif type(other) == type(0.0):
            z.handle = pygadgetron.cGT_axpby\
                (other, 0, self.handle, 0, 0, self.handle)
            return z;
        else:
            raise error('wrong multiplier')
    @staticmethod
    def axpby(a, x, b, y):
        '''
        Returns a linear combination a*x + b*y of two containers x and y.
        a and b: complex scalars
        x and y: DataContainers
        '''
        z = self.same_object()
        z.handle = pygadgetron.cGT_axpby\
            (a.real, a.imag, x.handle, b.real, b.imag, y.handle)
        return z;

class CoilImages(DataContainer):
    '''
    Class for a coil images container.
    Each item in the container is a 4D complex array of coil images values 
    on an xyz-slice (z-dimension is normally 1).
    '''
    def __init__(self):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject('CoilImagesList')
        check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteObject(self.handle)
    def calculate(self, acqs):
        '''
        Calculates coil images from a given sorted acquisitions.
        acqs: AcquisitionData
        '''
        if acqs.is_sorted() is False:
            print('WARNING: acquisitions may be in a wrong order')
        handle = pygadgetron.cGT_computeCoilImages\
            (self.handle, acqs.handle)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
    def image_dimensions(self):
        '''
        Returns each coil images array dimensions as a tuple (nx, ny, nz, nc),
        where nc is the number of active coils and nx, ny, nz are slice
        dimensions.
        '''
        dim = numpy.ndarray((4,), dtype = numpy.int32)
        pygadgetron.cGT_getCoilDataDimensions\
            (self.handle, 0, dim.ctypes.data)
        return tuple(numpy.asarray(dim))
    def as_array(self, ci_num):
        '''
        Returns specified coil images array as Numpy ndarray.
        ci_num: coil images array (slice) number
        '''
        nx, ny, nz, nc = self.image_dimensions()
        if nx == 0 or ny == 0 or nz == 0 or nc == 0:
            raise error('image data not available')
        re = numpy.ndarray((nc, nz, ny, nx), dtype = numpy.float64)
        im = numpy.ndarray((nc, nz, ny, nx), dtype = numpy.float64)
        pygadgetron.cGT_getCoilData\
            (self.handle, ci_num, re.ctypes.data, im.ctypes.data)
        return re + 1j * im

class CoilSensitivityMaps(DataContainer):
    '''
    Class for a coil sensitivity maps (csm) container.
    Each item in the container is a 4D complex array of csm values on an 
    xyz-slice (z-dimension is normally 1).
    '''
    def __init__(self):
        self.handle = None
        self.smoothness = 0
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def read(self, file):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
        self.handle = pygadgetron.cGT_CoilSensitivities(file)
        check_status(self.handle)
    def calculate(self, data, method = None):
        '''
        Calculates coil sensitivity maps from coil images or sorted 
        acquisitions.
        data  : either AcquisitionData or CoilImages
        method: either SRSS (Square Root of the Sum of Squares, default) or 
                Inati
        '''
        if isinstance(data, AcquisitionData):
            if data.is_sorted() is False:
                print('WARNING: acquisitions may be in a wrong order')
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
        self.handle = pygadgetron.cGT_CoilSensitivities('')
        check_status(self.handle)
        if method is not None:
            method_name, parm_list = name_and_parameters(method)
            parm = parse_arglist(parm_list)
        else:
            method_name = 'SRSS'
            parm = {}
        if isinstance(data, AcquisitionData):
            handle = pygadgetron.cGT_computeCoilSensitivities\
                (self.handle, data.handle)
            check_status(handle)
            pyiutil.deleteDataHandle(handle)
        elif isinstance(data, CoilImages):
            if method_name == 'Inati':
                if not HAVE_ISMRMRDTOOLS:
                    raise error('Inati method requires ismrmrd-python-tools')
                nz = data.number()
                for z in range(nz):
                    ci = numpy.squeeze(data.as_array(z))
                    (csm, rho) = coils.calculate_csm_inati_iter(ci)
                    self.append(csm)
            elif method_name == 'SRSS':
                if 'niter' in parm:
                    nit = int(parm['niter'])
                    _set_int_par\
                        (self.handle, 'coil_sensitivity', 'smoothness', nit)
                handle = pygadgetron.cGT_computeCSMsFromCIs\
                    (self.handle, data.handle)
                check_status(handle)
                pyiutil.deleteDataHandle(handle)
            else:
                raise error('Unknown method %s' % method_name)
    def append(self, csm):
        '''
        Appends a coil sensitivity map to self.
        csm: Numpy ndarray with csm values
        '''
        if self.handle is None:
            self.handle = pygadgetron.cGT_CoilSensitivities('')
            check_status(self.handle)
        shape = csm.shape
        nc = shape[0]
        if csm.ndim == 4:
            nz = shape[1]
            iy = 2
        else:
            nz = 1
            iy = 1
        ny = shape[iy]
        nx = shape[iy + 1]
        re = csm.real.copy()
        im = csm.imag.copy()
        handle = pygadgetron.cGT_appendCSM\
            (self.handle, nx, ny, nz, nc, re.ctypes.data, im.ctypes.data)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
    def map_dimensions(self):
        '''
        Returns each csm dimensions as a tuple (nx, ny, nz, nc),
        where nc is the number of active coils and nx, ny, nz are slice
        dimensions.
        '''
        dim = numpy.ndarray((4,), dtype = numpy.int32)
        pygadgetron.cGT_getCoilDataDimensions\
            (self.handle, 0, dim.ctypes.data)
        return tuple(numpy.asarray(dim))
    def as_array(self, csm_num):
        '''
        Returns specified csm as Numpy ndarray.
        csm_num: csm (slice) number
        '''
        nx, ny, nz, nc = self.map_dimensions()
        if nx == 0 or ny == 0 or nz == 0 or nc == 0:
            raise error('image data not available')
        re = numpy.ndarray((nc, nz, ny, nx), dtype = numpy.float64)
        im = numpy.ndarray((nc, nz, ny, nx), dtype = numpy.float64)
        pygadgetron.cGT_getCoilData\
            (self.handle, csm_num, re.ctypes.data, im.ctypes.data)
        return re + 1j * im
    def abs_as_array(self, csm_num):
        '''
        Returns the abs of specified csm as Numpy ndarray.
        csm_num: csm (slice) number
        '''
        nx, ny, nz, nc = self.map_dimensions()
        if nx == 0 or ny == 0 or nz == 0 or nc == 0:
            raise error('image data not available')
        array = numpy.ndarray((nc, nz, ny, nx), dtype = numpy.float64)
        pygadgetron.cGT_getCoilDataAbs\
            (self.handle, csm_num, array.ctypes.data)
        return array

class ImageData(DataContainer):
    '''
    Class for an MR images container.
    Each item in the container is a 3D complex or float array of the image 
    values on an xyz-slice (z-dimension is normally 1).
    '''
    def __init__(self):
        self.handle = None
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteObject(self.handle)
    def data_type(self, im_num):
        '''
        Returns the data type for a specified image (see 8 data types above).
        im_num: image (slice) 
        '''
        handle = pygadgetron.cGT_imageDataType(self.handle, im_num)
        check_status(handle)
        n = pyiutil.intDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return n
    def is_real(self):
        t = self.data_type(0)
        return t is not ISMRMRD_CXFLOAT and t is not ISMRMRD_CXDOUBLE
    def process(self, list):
        '''
        Processes self with an image processor specified by a list of gadgets.
        list: Python list of gadget description strings, each gadget 
              description being a string of the form
                '[label:]gadget_name[(property1=value1[,...])]'
              (square brackets embrace optional items, ... stands for etc.)
        '''
        ip = ImagesProcessor(list)
        return ip.process(self)
    def show(self):
        '''
        Interactively displays self's images.
        '''
        if not HAVE_PYLAB:
            print('pylab not found')
            return
        ni = self.number()
        if ni == 1:
            print('%d image' % ni)
        else:
            print('%d images' % ni)
        if ni < 1:
            return
        data = self.as_array()
        if not self.is_real():
            data = abs(data)
        print('Please enter the number of the image to view')
        print('(a value outside the range [1 : %d] will stop this loop)' % ni)
        while True:
            s = str(input('image: '))
            if len(s) < 1:
                break
            i = int(s)
            if i < 1 or i > ni:
                break
            pylab.figure(i)
            pylab.title('image %d' % i)
            pylab.imshow(data[i - 1, :, :])
            print('Close Figure %d window to continue...' % i)
            pylab.show()
    def write(self, out_file, out_group):
        '''
        Writes self's images to an hdf5 file.
        out_file : the file name (Python string)
        out_group: hdf5 dataset name (Python string)
        '''
        handle = pygadgetron.cGT_writeImages\
            (self.handle, out_file, out_group)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
    def select(self, attr, value):
        '''
        Creates an images container with images from self with the specified
        value of specified attribute.
        attr : the name of the attribute (Python string)
        value: the value of the attribute (Python string)
        '''
        images = ImageData()
        images.handle = pygadgetron.cGT_selectImages(self.handle, attr, value)
        check_status(images.handle)
        return images
    def as_array(self):
        '''
        Returns all self's images as a 3D Numpy ndarray.
        '''
        if self.number() < 1:
            return numpy.ndarray((0,0,0), dtype = numpy.float64)
        dim = numpy.ndarray((4,), dtype = numpy.int32)
        pygadgetron.cGT_getImageDimensions\
            (self.handle, 0, dim.ctypes.data)
        nx = dim[0]
        ny = dim[1]
        nz = dim[2]
        nc = dim[3]
        nz = nz*nc*self.number()
        if self.is_real():
            array = numpy.ndarray((nz, ny, nx), dtype = numpy.float64)
            pygadgetron.cGT_getImagesDataAsDoubleArray\
                (self.handle, array.ctypes.data)
            return array
        else:
            re = numpy.ndarray((nz, ny, nx), dtype = numpy.float64)
            im = numpy.ndarray((nz, ny, nx), dtype = numpy.float64)
            pygadgetron.cGT_getImagesDataAsComplexArray\
                (self.handle, re.ctypes.data, im.ctypes.data)
            return re + 1j*im

class AcquisitionInfo(PyGadgetronObject):
    '''
    Class for acquisition information parameters.
    '''
    def __init__(self):
        self.flags = 0
        self.encode_step_1 = 0
        self.slice = 0
        self.repetition = 0

class Acquisition(PyGadgetronObject):
    def __init__(self, file = None):
        self.handle = None
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteObject(self.handle)
    def flags(self):
        '''
        Returns acquisition flags as an integer (each bit corresponding to a 
        flag).
        '''
        return _int_par(self.handle, 'acquisition', 'flags')
    def number_of_samples(self):
        '''
        returns the number of samples in the readout direction.
        '''
        return _int_par(self.handle, 'acquisition', 'number_of_samples')
    def active_channels(self):
        '''
        Returns the number of active channels (coils).
        '''
        return _int_par(self.handle, 'acquisition', 'active_channels')
    def trajectory_dimensions(self):
        return _int_par(self.handle, 'acquisition', 'trajectory_dimensions')
    def idx_kspace_encode_step_1(self):
        return _int_par(self.handle, 'acquisition', 'idx_kspace_encode_step_1')
    def idx_repetition(self):
        return _int_par(self.handle, 'acquisition', 'idx_repetition')
    def idx_slice(self):
        return _int_par(self.handle, 'acquisition', 'idx_slice')

class AcquisitionData(DataContainer):
    '''
    Class for an MR acquisitions container.
    Each item is a 2D complex array of acquisition samples for each coil.
    '''
    def __init__(self, file = None):
        self.handle = None
        self.sorted = False
        self.info = None
        if file is not None:
            self.handle = pygadgetron.cGT_ISMRMRDAcquisitionsFromFile(file)
            check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteObject(self.handle)
    def number_of_acquisitions(self, select = 'all'):
        dim = self.dimensions(select)
        return dim[0]
    def sort(self):
        '''
        Sorts acquisitions with respect to (in this order):
            - repetition
            - slice
            - kspace_encode_step_1
        '''
        handle = pygadgetron.cGT_orderAcquisitions(self.handle)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
        self.sorted = True
    def is_sorted(self):
        return self.sorted
    def is_undersampled(self):
        return _int_par(self.handle, 'acquisitions', 'undersampled')
    def process(self, list):
        '''
        Processes self with an acquisition processor specified by a list of 
        gadgets.
        list: Python list of gadget description strings, each gadget 
              description being a string of the form
                '[label:]gadget_name[(property1=value1[,...])]'
              (square brackets embrace optional items, ... stands for etc.)
        '''
        ap = AcquisitionsProcessor(list)
        return ap.process(self)
    def acquisition(self, num):
        '''
        Returns the specified acquisition.
        num: acquisition number
        '''
        acq = Acquisition()
        acq.handle = pygadgetron.cGT_acquisitionFromContainer(self.handle, num)
        return acq
    def dimensions(self, select = 'all'):
        '''
        Returns acquisitions dimensions as a tuple (na, nc, ns), where na is
        the number of acquisitions, nc the number of coils and ns the number of
        samples.
        '''
        dim = numpy.ones((MAX_ACQ_DIMENSIONS,), dtype = numpy.int32)
        hv = pygadgetron.cGT_getAcquisitionsDimensions\
             (self.handle, dim.ctypes.data)
        #nr = pyiutil.intDataFromHandle(hv)
        pyiutil.deleteDataHandle(hv)
        if select == 'all':
            dim[2] = self.number()
        else:
            dim[2] = numpy.prod(dim[2:])
        return tuple(dim[2::-1])
    def set_info(self):
        '''
        Fills the array self.info with information for each acquisition.
        '''
        na, nc, ns = self.dimensions()
        self.info = numpy.empty((na,), dtype = object)
        for a in range(na):
            acq = self.acquisition(a)
            info = AcquisitionInfo()
            info.flags = acq.flags()
            info.encode_step_1 = acq.idx_kspace_encode_step_1()
            info.slice = acq.idx_slice()
            info.repetition = acq.idx_repetition()
            self.info[a] = info
    def get_info(self, par):
        '''
        Returns the array of values of the specified acquisition information 
        parameter.
        par: parameter name
        '''
        na, nc, ns = self.dimensions()
        if self.info is None:
            self.set_info()
        if par == 'flags':
            flags = numpy.empty((na,), dtype = numpy.int64)
            for a in range(na):
                flags[a] = self.info[a].flags
            return flags
        elif par == 'encode_step_1':
            es1 = numpy.empty((na,), dtype = numpy.int32)
            for a in range(na):
                es1[a] = self.info[a].encode_step_1
            return es1
        elif par == 'slice':
            s = numpy.empty((na,), dtype = numpy.int32)
            for a in range(na):
                s[a] = self.info[a].slice
            return s
        elif par == 'repetition':
            r = numpy.empty((na,), dtype = numpy.int32)
            for a in range(na):
                r[a] = self.info[a].repetition
            return r
        else:
            raise error('unknown acquisition parameter ' + par)
    def as_array(self, select = 'all'):
        '''
        Returns all self's acquisitions as a 3D Numpy ndarray.
        '''
        na = self.number()
        ny, nc, ns = self.dimensions(select)
        if select == 'all':
            n = na
        else:
            n = na + 1
        re = numpy.ndarray((ny, nc, ns), dtype = numpy.float64)
        im = numpy.ndarray((ny, nc, ns), dtype = numpy.float64)
        hv = pygadgetron.cGT_getAcquisitionsData\
            (self.handle, n, re.ctypes.data, im.ctypes.data)
        pyiutil.deleteDataHandle(hv)
        return re + 1j*im
    def fill(self, data):
        '''
        Fills self's acquisitions with specified values.
        data: Python Numpy array
        '''
        if self.handle is None:
            raise error('Undefined AcquisitionData object cannot be filled')
        na, nc, ns = data.shape
        re = numpy.copy(numpy.real(data))
        im = numpy.copy(numpy.imag(data))
        handle = pygadgetron.cGT_setAcquisitionsData\
            (self.handle, na, nc, ns, re.ctypes.data, im.ctypes.data)
        check_status(handle)
        pyiutil.deleteObject(self.handle)
        self.handle = handle

class AcquisitionModel(PyGadgetronObject):
    '''
    Class for MR acquisition model, an operator that maps images into
    simulated acquisitions.
    '''
    def __init__(self, acqs, imgs):
        self.handle = None
        self.handle = \
            pygadgetron.cGT_AcquisitionModel(acqs.handle, imgs.handle)
        check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteObject(self.handle)
    def set_coil_sensitivity_maps(self, csm):
        '''
        Specifies the coil sensitivity maps to be used by the model.
        csm: CoilSensitivityMaps
        '''
        handle = pygadgetron.cGT_setCSMs(self.handle, csm.handle)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
        self.csms = 'set'
    def forward(self, images):
        acqs = AcquisitionData()
        acqs.handle = pygadgetron.cGT_AcquisitionModelForward\
            (self.handle, images.handle)
        check_status(acqs.handle)
        return acqs;
    def backward(self, acqs):
        images = ImageData()
        images.handle = pygadgetron.cGT_AcquisitionModelBackward\
            (self.handle, acqs.handle)
        check_status(images.handle)
        return images

class Gadget(PyGadgetronObject):
    def __init__(self, name):
        self.handle = None
        name, prop = name_and_parameters(name)
        self.handle = pygadgetron.cGT_newObject(name)
        check_status(self.handle)
        if prop is not None:
            self.set_properties(prop)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteObject(self.handle)
    def set_property(self, prop, value):
        handle = pygadgetron.cGT_setGadgetProperty(self.handle, prop, value)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
    def set_properties(self, prop):
        handle = pygadgetron.cGT_setGadgetProperties(self.handle, prop)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
    def value_of(self, prop):
        return _char_par(self.handle, 'gadget', prop)

class GadgetChain(PyGadgetronObject):
    def __init__(self):
        self.handle = pygadgetron.cGT_newObject('GadgetChain')
        check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteObject(self.handle)
    def add_reader(self, id, reader):
        handle = pygadgetron.cGT_addReader(self.handle, id, reader.handle)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
    def add_writer(self, id, writer):
        handle = pygadgetron.cGT_addWriter(self.handle, id, writer.handle)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
    def add_gadget(self, id, gadget):
        handle = pygadgetron.cGT_addGadget(self.handle, id, gadget.handle)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
    def set_gadget_property(self, id, prop, value):
        if type(value) == type('abc'):
            v = value
        else:
            v = repr(value).lower()
        hg = _parameterHandle(self.handle, 'gadget_chain', id)
        handle = pygadgetron.cGT_setGadgetProperty(hg, prop, v)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
        pyiutil.deleteDataHandle(hg)
    def value_of_gadget_property(self, id, prop):
        hg = _parameterHandle(self.handle, 'gadget_chain', id)
        hv = _parameterHandle(hg, 'gadget', prop)
        value = pyiutil.charDataFromHandle(hv)
        pyiutil.deleteDataHandle(hg)
        pyiutil.deleteDataHandle(hv)
        return value

class ImagesReconstructor(GadgetChain):
    def __init__(self, list = None):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject('ImagesReconstructor')
        check_status(self.handle)
        self.input_data = None
        if list is None:
            return
        for i in range(len(list)):
            label, name = label_and_name(list[i])
            self.add_gadget(label, Gadget(name))
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def set_input(self, input_data):
        self.input_data = input_data
    def process(self):
        if self.input_data is None:
            raise error('no input data')
        handle = pygadgetron.cGT_reconstructImages\
             (self.handle, self.input_data.handle)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
    def get_output(self, subset = None):
        output = ImageData()
        output.handle = pygadgetron.cGT_reconstructedImages(self.handle)
        check_status(output.handle)
        if subset is None:
            return output
        else:
            return output.select('GADGETRON_DataRole', subset)
    def reconstruct(self, input_data):
        handle = pygadgetron.cGT_reconstructImages\
             (self.handle, input_data.handle)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
        images = ImageData()
        images.handle = pygadgetron.cGT_reconstructedImages(self.handle)
        check_status(images.handle)
        return images

class ImagesProcessor(GadgetChain):
    def __init__(self, list = None):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject('ImagesProcessor')
        check_status(self.handle)
        self.input_data = None
        if list is None:
            return
        for i in range(len(list)):
            self.add_gadget('g' + repr(i + 1), Gadget(list[i]))
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteObject(self.handle)
##    def set_input(self, input_data):
##        self.input_data = input_data
    def process(self, input_data):
##        if self.input_data is None:
##            raise error('no input data')
        images = ImageData()
        images.handle = pygadgetron.cGT_processImages\
             (self.handle, input_data.handle)
        check_status(images.handle)
        return images
##    def get_output(self):
##        images = ImageData()
##        images.handle = pygadgetron.cGT_reconstructedImagesList(self.handle)
##        check_status(images.handle)
##        return images

class AcquisitionsProcessor(GadgetChain):
    def __init__(self, list = None):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject('AcquisitionsProcessor')
        check_status(self.handle)
        self.input_data = None
        if list is None:
            return
        for i in range(len(list)):
            self.add_gadget('g' + repr(i + 1), Gadget(list[i]))
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteObject(self.handle)
    def process(self, input_data):
        acquisitions = AcquisitionData()
        acquisitions.handle = pygadgetron.cGT_processAcquisitions\
             (self.handle, input_data.handle)
        check_status(acquisitions.handle)
        return acquisitions

class SimpleReconstruction(ImagesReconstructor):
    def __init__(self):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject('SimpleReconstructionProcessor')
        check_status(self.handle)
        self.input_data = None
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteObject(self.handle)
    
class GenericCartesianGRAPPAReconstruction(ImagesReconstructor):
    def __init__(self):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject\
            ('SimpleGRAPPAReconstructionProcessor')
        check_status(self.handle)
        self.input_data = None
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteObject(self.handle)
    def compute_gfactors(self, flag):
        self.set_gadget_property('gadget4', 'send_out_gfactor', flag)
    
def MR_remove_x_oversampling(input_data):
    handle = pygadgetron.cGT_newObject('RemoveOversamplingProcessor')
    check_status(handle)
    output_data = AcquisitionData()
    output_data.handle = pygadgetron.cGT_processAcquisitions\
         (handle, input_data.handle)
    check_status(output_data.handle)
    pyiutil.deleteObject(handle)
    return output_data

def preprocess_acquisitions(input_data):
    return input_data.process(\
        ['NoiseAdjustGadget', \
         'AsymmetricEchoAdjustROGadget', \
         'RemoveROOversamplingGadget'])
    
