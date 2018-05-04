''' 
Object-Oriented wrap for the cGadgetron-to-Python interface pygadgetron.py
'''

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
##
## This is software developed for the Collaborative Computational
## Project in Positron Emission Tomography and Magnetic Resonance imaging
## (http://www.ccppetmr.ac.uk/).
##
## Licensed under the Apache License, Version 2.0 (the "License");
##   you may not use this file except in compliance with the License.
##   You may obtain a copy of the License at
##       http://www.apache.org/licenses/LICENSE-2.0
##   Unless required by applicable law or agreed to in writing, software
##   distributed under the License is distributed on an "AS IS" BASIS,
##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
##   See the License for the specific language governing permissions and
##   limitations under the License.

import abc
import numpy
import os
try:
    import pylab
    HAVE_PYLAB = True
except:
    HAVE_PYLAB = False
import sys
import time
##try:
##    from ismrmrdtools import coils
##    HAVE_ISMRMRDTOOLS = True
##except:
##    print('ismrmrd-python-tools not installed')
##    HAVE_ISMRMRDTOOLS = False

from pUtilities import *
import pyiutilities as pyiutil
import pygadgetron

if sys.version_info[0] >= 3 and sys.version_info[1] >= 4:
    ABC = abc.ABC
else:
    ABC = abc.ABCMeta('ABC', (), {})

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
    try_calling(pygadgetron.cGT_setParameter(hs, set, par, hv))
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
def _int_pars(handle, set, par, n):
    h = pygadgetron.cGT_parameter(handle, set, par)
    check_status(h)
    for i in range(n):
        value += (pyiutil.intDataItemFromHandle(h, n - 1 - i),)
    pyiutil.deleteDataHandle(h)
    return value
def _uint16_pars(handle, set, par, n):
    h = pygadgetron.cGT_parameter(handle, set, par)
    check_status(h)
    #value = (pyiutil.uint16DataItemFromHandle(h, 0),)
    #for i in range(1, n):
    value = ()
    for i in range(n):
        value += (pyiutil.uint16DataItemFromHandle(h, n - 1 - i),)
    pyiutil.deleteDataHandle(h)
    return value
def _uint32_pars(handle, set, par, n):
    h = pygadgetron.cGT_parameter(handle, set, par)
    check_status(h)
    value = ()
    for i in range(n):
        value += (pyiutil.uint32DataItemFromHandle(h, n - 1 - i),)
    pyiutil.deleteDataHandle(h)
    return value
def _uint64_pars(handle, set, par, n):
    h = pygadgetron.cGT_parameter(handle, set, par)
    check_status(h)
    value = ()
    for i in range(n):
        value += (pyiutil.uint64DataItemFromHandle(h, n - 1 - i),)
    pyiutil.deleteDataHandle(h)
    return value
def _char_par(handle, set, par):
    h = pygadgetron.cGT_parameter(handle, set, par)
    check_status(h)
    value = pyiutil.charDataFromHandle(h)
    pyiutil.deleteDataHandle(h)
    return value
def _float_par(handle, set, par):
    h = pygadgetron.cGT_parameter(handle, set, par)
    check_status(h)
    v = pyiutil.floatDataFromHandle(h)
    pyiutil.deleteDataHandle(h)
    return v
def _float_pars(handle, set, par, n):
    h = pygadgetron.cGT_parameter(handle, set, par)
    check_status(h)
    value = ()
    for i in range(n):
        value += (pyiutil.floatDataItemFromHandle(h, n - 1 - i),)
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

### low-level client functionality
### likely to be obsolete- not used for a long time
##class ClientConnector:
##    def __init__(self):
##        self.handle = None
##        self.handle = pygadgetron.cGT_newObject('GTConnector')
##        check_status(self.handle)
##    def __del__(self):
##        if self.handle is not None:
##            pyiutil.deleteDataHandle(self.handle)
##    def set_timeout(self, timeout):
##        handle = pygadgetron.cGT_setConnectionTimeout(self.handle, timeout)
##        check_status(handle)
##        pyiutil.deleteDataHandle(handle)
##    def connect(self, host, port):
##        handle = pygadgetron.cGT_connect(self.handle, host, port)
##        check_status(handle)
##        pyiutil.deleteDataHandle(handle)
##    def disconnect(self):
##        handle = pygadgetron.cGT_disconnect(self.handle)
##        check_status(handle)
##        pyiutil.deleteDataHandle(handle)
##    def register_images_receiver(self, imgs):
##        handle = pygadgetron.cGT_registerImagesReceiver\
##            (self.handle, imgs.handle)
##        check_status(handle)
##        pyiutil.deleteDataHandle(handle)
##    def config_gadget_chain(self, gc):
##        handle = pygadgetron.cGT_configGadgetChain(self.handle, gc.handle)
##        check_status(handle)
##        pyiutil.deleteDataHandle(handle)
##    def send_config_file(self, file):
##        handle = pygadgetron.cGT_sendConfigFile(self.handle, file)
##        check_status(handle)
##        pyiutil.deleteDataHandle(handle)
##    def send_parameters(self, par):
##        handle = pygadgetron.cGT_sendParameters(self.handle, par)
##        check_status(handle)
##        pyiutil.deleteDataHandle(handle)
##    def send_acquisitions(self, acq):
##        handle = pygadgetron.cGT_sendAcquisitions(self.handle, acq.handle)
##        check_status(handle)
##        pyiutil.deleteDataHandle(handle)
##    def send_images(self, img):
##        handle = pygadgetron.cGT_sendImages(self.handle, img.handle)
##        check_status(handle)
##        pyiutil.deleteDataHandle(handle)

# base class for complex data container classes
class DataContainer(ABC):
    '''
    Class for an abstract data container.
    '''
    def __init__(self):
        self.handle = None
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    @abc.abstractmethod
    def same_object(self):
        '''
        Returns an object of the same type as self.
        '''
        pass
    def number(self):
        '''
        Returns the number of items in the container.
        '''
        assert self.handle is not None
        handle = pygadgetron.cGT_dataItems(self.handle)
        check_status(handle)
        n = pyiutil.intDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return n
    def norm(self):
        '''
        Returns the 2-norm of the container data viewed as a vector.
        '''
        assert self.handle is not None
        handle = pygadgetron.cGT_norm(self.handle)
        check_status(handle)
        r = pyiutil.floatDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return r;
    def dot(self, other):
        '''
        Returns the dot product of the container data with another container 
        data viewed as vectors.
        other: DataContainer
        '''
##        assert self.handle is not None
##        assert type(self) == type(other)
        assert_validities(self, other)
        handle = pygadgetron.cGT_dot(self.handle, other.handle)
        check_status(handle)
        re = pyiutil.floatReDataFromHandle(handle)
        im = pyiutil.floatImDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return complex(re, im)
    def __add__(self, other):
        '''
        Overloads + for data containers.
        Returns the sum of the container data with another container 
        data viewed as vectors.
        other: DataContainer
        '''
##        assert self.handle is not None
##        assert type(self) == type(other)
        assert_validities(self, other)
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
##        assert self.handle is not None
##        assert type(self) == type(other)
        assert_validities(self, other)
        z = self.same_object()
        z.handle = pygadgetron.cGT_axpby\
            (1.0, 0.0, self.handle, -1.0, 0.0, other.handle)
        return z;
    def __mul__(self, other):
        '''
        Overloads * for data containers multiplication by a scalar or another
        data container. Returns the product self*other if other is a scalar
        or the dot product if it is DataContainer.
        other: DataContainer or a (real or complex) scalar
        '''
        assert self.handle is not None
        if type(self) == type(other):
            return self.dot(other)
        z = self.same_object()
        if type(other) == type(0):
            other = float(other)
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
        Overloads * for data containers multiplication by a scalar from
        the left, i.e. computes and returns the product other*self.
        other: a real or complex scalar
        '''
        assert self.handle is not None
        z = self.same_object()
        if type(other) == type(0):
            other = float(other)
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
##        assert type(x) == type(y)
        assert_validities(x, y)
        z = x.same_object()
        z.handle = pygadgetron.cGT_axpby\
            (a.real, a.imag, x.handle, b.real, b.imag, y.handle)
        return z;

class CoilImageData(DataContainer):
    '''
    Class for a coil images container.
    Each item in the container is a 4D complex array of coil images values 
    on an xyz-slice (z-dimension is normally 1).
    '''
    def __init__(self):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject('CoilImages')
        check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def same_object(self):
        return CoilImageData()
    def calculate(self, acqs):
        '''
        Calculates coil images from a given sorted acquisitions.
        acqs: AcquisitionData
        '''
        #assert isinstance(acqs, AcquisitionData)
        assert_validity(acqs, AcquisitionData)
        if acqs.is_sorted() is False:
            print('WARNING: acquisitions may be in a wrong order')
        try_calling(pygadgetron.cGT_computeCoilImages\
            (self.handle, acqs.handle))
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
        re = numpy.ndarray((nc, nz, ny, nx), dtype = numpy.float32)
        im = numpy.ndarray((nc, nz, ny, nx), dtype = numpy.float32)
        pygadgetron.cGT_getCoilData\
            (self.handle, ci_num, re.ctypes.data, im.ctypes.data)
        return re + 1j * im

DataContainer.register(CoilImageData)

class CoilSensitivityData(DataContainer):
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
    def same_object(self):
        return CoilSensitivityData()
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
            assert data.handle is not None
            _set_int_par\
                (self.handle, 'coil_sensitivity', 'smoothness', self.smoothness)
            try_calling(pygadgetron.cGT_computeCoilSensitivities\
                (self.handle, data.handle))
        elif isinstance(data, CoilImageData):
            assert data.handle is not None
            if method_name == 'Inati':
#                if not HAVE_ISMRMRDTOOLS:
                try:
                    from ismrmrdtools import coils
                except:
                    raise error('Inati method requires ismrmrd-python-tools')
                nz = data.number()
                for z in range(nz):
                    ci = numpy.squeeze(data.as_array(z))
                    (csm, rho) = coils.calculate_csm_inati_iter(ci)
                    self.append(csm.astype(numpy.complex64))
            elif method_name == 'SRSS':
                if 'niter' in parm:
                    nit = int(parm['niter'])
                    _set_int_par\
                        (self.handle, 'coil_sensitivity', 'smoothness', nit)
                try_calling(pygadgetron.cGT_computeCSMsFromCIs\
                    (self.handle, data.handle))
            else:
                raise error('Unknown method %s' % method_name)
        else:
            raise error('Cannot calculate coil sensitivities from %s' % \
                        repr(type(data)))
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
        assert self.handle is not None
        dim = numpy.ndarray((4,), dtype = numpy.int32)
        pygadgetron.cGT_getCoilDataDimensions\
            (self.handle, 0, dim.ctypes.data)
        return tuple(numpy.asarray(dim))
    def as_array(self, csm_num):
        '''
        Returns specified csm as Numpy ndarray.
        csm_num: csm (slice) number
        '''
        assert self.handle is not None
        nx, ny, nz, nc = self.map_dimensions()
        if nx == 0 or ny == 0 or nz == 0 or nc == 0:
            raise error('image data not available')
        re = numpy.ndarray((nc, nz, ny, nx), dtype = numpy.float32)
        im = numpy.ndarray((nc, nz, ny, nx), dtype = numpy.float32)
        pygadgetron.cGT_getCoilData\
            (self.handle, csm_num, re.ctypes.data, im.ctypes.data)
        return re + 1j * im
    def abs_as_array(self, csm_num):
        '''
        Returns the abs of specified csm as Numpy ndarray.
        csm_num: csm (slice) number
        '''
        assert self.handle is not None
        nx, ny, nz, nc = self.map_dimensions()
        if nx == 0 or ny == 0 or nz == 0 or nc == 0:
            raise error('image data not available')
        array = numpy.ndarray((nc, nz, ny, nx), dtype = numpy.float32)
        pygadgetron.cGT_getCoilDataAbs\
            (self.handle, csm_num, array.ctypes.data)
        return array

DataContainer.register(CoilSensitivityData)

class Image:
    '''
    Class for an MR image.
    '''
    def __init__(self, image_data = None, image_num = 0):
        self.handle = None
        if image_data is not None:
            self.handle = pygadgetron.cGT_imageWrapFromContainer \
                          (image_data.handle, image_num)
            check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def is_real(self):
        assert self.handle is not None
        t = self.data_type()
        return t is not ISMRMRD_CXFLOAT and t is not ISMRMRD_CXDOUBLE
    def as_array(self):
        '''
        Returns image data as a 3D Numpy ndarray.
        '''
        assert self.handle is not None
        dim = numpy.ndarray((4,), dtype = numpy.int32)
        pygadgetron.cGT_getImageDim(self.handle, dim.ctypes.data)
        nx = dim[0]
        ny = dim[1]
        nz = dim[2]
        nc = dim[3]
        nz = nz*nc
        if self.is_real():
            array = numpy.ndarray((nz, ny, nx), dtype = numpy.float32)
            pygadgetron.cGT_getImageDataAsFloatArray\
                (self.handle, array.ctypes.data)
            return array
        else:
            re = numpy.ndarray((nz, ny, nx), dtype = numpy.float32)
            im = numpy.ndarray((nz, ny, nx), dtype = numpy.float32)
            pygadgetron.cGT_getImageDataAsComplexArray\
                (self.handle, re.ctypes.data, im.ctypes.data)
            return re + 1j*im
    def version(self):
        assert self.handle is not None
        return _int_par(self.handle, 'image', 'version')
    def flags(self):
        assert self.handle is not None
        return _int_par(self.handle, 'image', 'flags')
    def data_type(self):
        assert self.handle is not None
        return _int_par(self.handle, 'image', 'data_type')
    def measurement_uid(self):
        assert self.handle is not None
        return _int_par(self.handle, 'image', 'measurement_uid')
    def channels(self):
        assert self.handle is not None
        return _int_par(self.handle, 'image', 'channels')
    def average(self):
        assert self.handle is not None
        return _int_par(self.handle, 'image', 'average')
    def slice(self):
        assert self.handle is not None
        return _int_par(self.handle, 'image', 'slice')
    def contrast(self):
        assert self.handle is not None
        return _int_par(self.handle, 'image', 'contrast')
    def phase(self):
        assert self.handle is not None
        return _int_par(self.handle, 'image', 'phase')
    def repetition(self):
        assert self.handle is not None
        return _int_par(self.handle, 'image', 'repetition')
    def set(self):
        assert self.handle is not None
        return _int_par(self.handle, 'image', 'set')
    def acquisition_time_stamp(self):
        assert self.handle is not None
        return _int_par(self.handle, 'image', 'acquisition_time_stamp')
    def image_type(self):
        assert self.handle is not None
        return _int_par(self.handle, 'image', 'image_type')
    def image_index(self):
        assert self.handle is not None
        return _int_par(self.handle, 'image', 'image_index')
    def image_series_index(self):
        assert self.handle is not None
        return _int_par(self.handle, 'image', 'image_series_index')
    def attribute_string_len(self):
        assert self.handle is not None
        return _int_par(self.handle, 'image', 'attribute_string_len')
    def matrix_size(self):
        assert self.handle is not None
        return _uint16_pars(self.handle, 'image', 'matrix_size', 3)
    def physiology_time_stamp(self):
        assert self.handle is not None
        return _uint32_pars(self.handle, 'image', 'physiology_time_stamp', 3)
    def field_of_view(self):
        assert self.handle is not None
        return _float_pars(self.handle, 'image', 'field_of_view', 3)
    def position(self):
        assert self.handle is not None
        return _float_pars(self.handle, 'image', 'position', 3)
    def read_dir(self):
        assert self.handle is not None
        return _float_pars(self.handle, 'image', 'read_dir', 3)
    def phase_dir(self):
        assert self.handle is not None
        return _float_pars(self.handle, 'image', 'phase_dir', 3)
    def slice_dir(self):
        assert self.handle is not None
        return _float_pars(self.handle, 'image', 'slice_dir', 3)
    def patient_table_position(self):
        assert self.handle is not None
        return _float_pars \
               (self.handle, 'image', 'patient_table_position', 3)
    def info(self, method):
        return eval('self.' + method + '()')

class ImageData(DataContainer):
    '''
    Class for an MR images container.
    Each item in the container is a 3D complex or float array of the image 
    values on an xyz-slice (z-dimension is normally 1).
    '''
    def __init__(self, file = None):
        self.handle = None
        if file is None:
            return
        self.handle = pygadgetron.cGT_readImages(file)
        check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def same_object(self):
        return ImageData()
    def read_from_file(self, file):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
        self.handle = pygadgetron.cGT_readImages(file)
        check_status(self.handle)
    def data_type(self, im_num):
        '''
        Returns the data type for a specified image (see 8 data types above).
        im_num: image (slice) 
        '''
        assert self.handle is not None
        handle = pygadgetron.cGT_imageDataType(self.handle, im_num)
        check_status(handle)
        n = pyiutil.intDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return n
    def is_real(self):
        assert self.handle is not None
        t = self.data_type(0)
        return t is not ISMRMRD_CXFLOAT and t is not ISMRMRD_CXDOUBLE
    def process(self, list):
        '''
        Returns processed self with an image processor specified by
        a list of gadgets.
        list: Python list of gadget description strings, each gadget 
              description being a string of the form
                '[label:]gadget_name[(property1=value1[,...])]'
              (square brackets embrace optional items, ... stands for etc.)
        '''
        assert self.handle is not None
        ip = ImageDataProcessor(list)
        return ip.process(self)
    def clone(self):
        '''
        Returns a copy of self.
        '''
        assert self.handle is not None
        ip = ImageDataProcessor()
        return ip.process(self)
    def image(self, im_num):
        return Image(self, im_num)
##    def show(self):
##        '''
##        Interactively displays self's images.
##        '''
##        assert self.handle is not None
##        if not HAVE_PYLAB:
##            print('pylab not found')
##            return
##        ni = self.number()
##        if ni == 1:
##            print('%d image' % ni)
##        else:
##            print('%d images' % ni)
##        if ni < 1:
##            return
##        data = self.as_array()
##        if not self.is_real():
##            data = abs(data)
##        print('Please enter the number of the image to view')
##        print('(a value outside the range [1 : %d] will stop this loop)' % ni)
##        while True:
##            s = str(input('image: '))
##            if len(s) < 1:
##                break
##            i = int(s)
##            if i < 1 or i > ni:
##                break
##            image = self.image(i - 1)
##            pylab.figure(i)
##            pylab.title('image %d' % i)
##            arr = image.as_array()
##            if not image.is_real():
##                arr = abs(arr)
##            pylab.imshow(arr[0,:,:])
##            #pylab.imshow(data[i - 1, :, :])
##            print('Close Figure %d window to continue...' % i)
##            pylab.show()
    def write(self, out_file, out_group):
        '''
        Writes self's images to an hdf5 file.
        out_file : the file name (Python string)
        out_group: hdf5 dataset name (Python string)
        '''
        assert self.handle is not None
        try_calling(pygadgetron.cGT_writeImages\
                    (self.handle, out_file, out_group))
    def select(self, attr, value):
        '''
        Creates an images container with images from self with the specified
        value of specified attribute.
        attr : the name of the attribute (Python string)
        value: the value of the attribute (Python string)
        '''
        assert self.handle is not None
        images = ImageData()
        images.handle = pygadgetron.cGT_selectImages(self.handle, attr, value)
        check_status(images.handle)
        return images
    def get_info(self, par):
        '''
        Returns the array of values of the specified image information 
        parameter.
        par: parameter name
        '''
        ni = self.number()
        info = numpy.empty((ni,), dtype = object)
        for i in range(ni):
            image = self.image(i)
            info[i] = image.info(par)
        return info
    def as_array(self):
        '''
        Returns all self's images as a 3D Numpy ndarray.
        '''
        assert self.handle is not None
        if self.number() < 1:
            return numpy.ndarray((0,0,0), dtype = numpy.float32)
        dim = numpy.ndarray((4,), dtype = numpy.int32)
        image = Image(self)
        pygadgetron.cGT_getImageDim(image.handle, dim.ctypes.data)
##        pygadgetron.cGT_getImageDimensions\
##            (self.handle, 0, dim.ctypes.data)
        nx = dim[0]
        ny = dim[1]
        nz = dim[2]
        nc = dim[3]
        nz = nz*nc*self.number()
        if self.is_real():
            array = numpy.ndarray((nz, ny, nx), dtype = numpy.float32)
            pygadgetron.cGT_getImagesDataAsFloatArray\
                (self.handle, array.ctypes.data)
            return array
        else:
            re = numpy.ndarray((nz, ny, nx), dtype = numpy.float32)
            im = numpy.ndarray((nz, ny, nx), dtype = numpy.float32)
            pygadgetron.cGT_getImagesDataAsComplexArray\
                (self.handle, re.ctypes.data, im.ctypes.data)
            return re + 1j*im
    def fill(self, data):
        '''
        Fills self's image data with specified values.
        data: Python Numpy array
        '''
        assert self.handle is not None
        re = numpy.real(data).astype(numpy.float32)
        im = numpy.imag(data).astype(numpy.float32)
        try_calling(pygadgetron.cGT_setComplexImagesData\
            (self.handle, re.ctypes.data, im.ctypes.data))

DataContainer.register(ImageData)

class Acquisition:
    def __init__(self, file = None):
        self.handle = None
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def version(self):
        assert self.handle is not None
        return _int_par(self.handle, 'acquisition', 'version')
    def flags(self):
        '''
        Returns acquisition flags as an integer (each bit corresponding to a 
        flag).
        '''
        assert self.handle is not None
        return _int_par(self.handle, 'acquisition', 'flags')
    def measurement_uid(self):
        assert self.handle is not None
        return _int_par(self.handle, 'acquisition', 'measurement_uid')
    def scan_counter(self):
        assert self.handle is not None
        return _int_par(self.handle, 'acquisition', 'scan_counter')
    def acquisition_time_stamp(self):
        assert self.handle is not None
        return _int_par(self.handle, 'acquisition', 'acquisition_time_stamp')
    def number_of_samples(self):
        '''
        returns the number of samples in the readout direction.
        '''
        assert self.handle is not None
        return _int_par(self.handle, 'acquisition', 'number_of_samples')
    def available_channels(self):
        assert self.handle is not None
        return _int_par(self.handle, 'acquisition', 'available_channels')
    def active_channels(self):
        '''
        Returns the number of active channels (coils).
        '''
        assert self.handle is not None
        return _int_par(self.handle, 'acquisition', 'active_channels')
    def discard_pre(self):
        assert self.handle is not None
        return _int_par(self.handle, 'acquisition', 'discard_pre')
    def discard_post(self):
        assert self.handle is not None
        return _int_par(self.handle, 'acquisition', 'discard_post')
    def center_sample(self):
        assert self.handle is not None
        return _int_par(self.handle, 'acquisition', 'center_sample')
    def encoding_space_ref(self):
        assert self.handle is not None
        return _int_par(self.handle, 'acquisition', 'encoding_space_ref')
    def trajectory_dimensions(self):
        assert self.handle is not None
        return _int_par(self.handle, 'acquisition', 'trajectory_dimensions')
    def kspace_encode_step_1(self):
        assert self.handle is not None
        return _int_par(self.handle, 'acquisition', 'idx_kspace_encode_step_1')
    def kspace_encode_step_2(self):
        assert self.handle is not None
        return _int_par(self.handle, 'acquisition', 'idx_kspace_encode_step_2')
    def average(self):
        assert self.handle is not None
        return _int_par(self.handle, 'acquisition', 'idx_average')
    def slice(self):
        assert self.handle is not None
        return _int_par(self.handle, 'acquisition', 'idx_slice')
    def contrast(self):
        assert self.handle is not None
        return _int_par(self.handle, 'acquisition', 'idx_contrast')
    def phase(self):
        assert self.handle is not None
        return _int_par(self.handle, 'acquisition', 'idx_phase')
    def repetition(self):
        assert self.handle is not None
        return _int_par(self.handle, 'acquisition', 'idx_repetition')
    def set(self):
        assert self.handle is not None
        return _int_par(self.handle, 'acquisition', 'idx_set')
    def segment(self):
        assert self.handle is not None
        return _int_par(self.handle, 'acquisition', 'idx_segment')
    def physiology_time_stamp(self):
        assert self.handle is not None
        return _uint32_pars(self.handle, 'acquisition', 'physiology_time_stamp', 3)
    def channel_mask(self):
        assert self.handle is not None
        return _uint64_pars(self.handle, 'acquisition', 'channel_mask', 16)
    def sample_time_us(self):
        assert self.handle is not None
        return _float_par(self.handle, 'acquisition', 'sample_time_us')
    def position(self):
        assert self.handle is not None
        return _float_pars(self.handle, 'acquisition', 'position', 3)
    def read_dir(self):
        assert self.handle is not None
        return _float_pars(self.handle, 'acquisition', 'read_dir', 3)
    def phase_dir(self):
        assert self.handle is not None
        return _float_pars(self.handle, 'acquisition', 'phase_dir', 3)
    def slice_dir(self):
        assert self.handle is not None
        return _float_pars(self.handle, 'acquisition', 'slice_dir', 3)
    def patient_table_position(self):
        assert self.handle is not None
        return _float_pars \
               (self.handle, 'acquisition', 'patient_table_position', 3)
    def info(self, method):
        return eval('self.' + method + '()')

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
            pyiutil.deleteDataHandle(self.handle)
    @staticmethod
    def set_storage_scheme(scheme):
        '''Sets acquisition data storage scheme.

        scheme = 'file' (default):
            all acquisition data generated from now on will be kept in
            scratch files deleted after the user's script terminates
        scheme = 'memory':
            all acquisition data generated from now on will be kept in RAM
            (avoid if data is very large)
        '''
        try_calling(pygadgetron.cGT_setAcquisitionsStorageScheme(scheme))
    @staticmethod
    def get_storage_scheme():
        '''Returns acquisition data storage scheme.
        '''
        handle = pygadgetron.cGT_getAcquisitionsStorageScheme()
        check_status(handle)
        scheme = pyiutil.charDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return scheme
    def same_object(self):
        return AcquisitionData()
##    def number_of_acquisitions(self, select = 'image'):
##        assert self.handle is not None
##        dim = self.dimensions(select)
##        return dim[0]
    def number_of_readouts(self, select = 'image'):
        dim = self.dimensions(select)
        return dim[0]
    def sort(self):
        '''
        Sorts acquisitions with respect to (in this order):
            - repetition
            - slice
            - kspace_encode_step_1
        '''
        assert self.handle is not None
        try_calling(pygadgetron.cGT_orderAcquisitions(self.handle))
        self.sorted = True
    def is_sorted(self):
        return self.sorted
    def is_undersampled(self):
        assert self.handle is not None
        return _int_par(self.handle, 'acquisitions', 'undersampled')
    def process(self, list):
        '''
        Returns processed self with an acquisition processor specified by
        a list of gadgets.
        list: Python list of gadget description strings, each gadget 
              description being a string of the form
                '[label:]gadget_name[(property1=value1[,...])]'
              (square brackets embrace optional items, ... stands for etc.)
        '''
        ap = AcquisitionDataProcessor(list)
        return ap.process(self)
    def clone(self):
        '''
        Returns a copy of self.
        '''
        ap = AcquisitionDataProcessor()
        return ap.process(self)
    def acquisition(self, num):
        '''
        Returns the specified acquisition.
        num: acquisition number
        '''
        assert self.handle is not None
        acq = Acquisition()
        acq.handle = pygadgetron.cGT_acquisitionFromContainer(self.handle, num)
        return acq
    def dimensions(self, select = 'image'):
        '''
        Returns acquisitions dimensions as a tuple (na, nc, ns), where na is
        the number of acquisitions, nc the number of coils and ns the number of
        samples.
        If select is set to 'all', the total number of acquisitions is returned.
        Otherwise, the number of acquisitions directly related to imaging data
        is returned.
        '''
        assert self.handle is not None
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
    def get_info(self, par):
        '''
        Returns the array of values of the specified acquisition information 
        parameter.
        par: parameter name
        '''
        na, nc, ns = self.dimensions()
        info = numpy.empty((na,), dtype = object)
        for a in range(na):
            acq = self.acquisition(a)
            info[a] = acq.info(par)
        return info
    def as_array(self, select = 'image'):
        '''
        Returns selected self's acquisitions as a 3D Numpy ndarray.
        '''
        assert self.handle is not None
        na = self.number()
        ny, nc, ns = self.dimensions(select)
        if select == 'all': # return all
            n = na
        else: # return only image-related
            n = na + 1
        re = numpy.ndarray((ny, nc, ns), dtype = numpy.float32)
        im = numpy.ndarray((ny, nc, ns), dtype = numpy.float32)
        hv = pygadgetron.cGT_getAcquisitionsData\
            (self.handle, n, re.ctypes.data, im.ctypes.data)
        pyiutil.deleteDataHandle(hv)
        return re + 1j*im
    def fill(self, data):
        '''
        Fills self's acquisitions with specified values.
        data: Python Numpy array
        '''
        assert self.handle is not None
        na, nc, ns = data.shape
        re = numpy.real(data).astype(numpy.float32)
        im = numpy.imag(data).astype(numpy.float32)
        try_calling(pygadgetron.cGT_setAcquisitionsData\
            (self.handle, na, nc, ns, re.ctypes.data, im.ctypes.data))
    def write(self, out_file):
        '''
        Writes self's acquisitions to an hdf5 file.
        out_file : the file name (Python string)
        '''
        assert self.handle is not None
        try_calling(pygadgetron.cGT_writeAcquisitions\
                    (self.handle, out_file))

DataContainer.register(AcquisitionData)

class AcquisitionModel:
    '''
    Class for MR acquisition model, an operator that maps images into
    simulated acquisitions.
    '''
    def __init__(self, acqs, imgs):
##        assert isinstance(acqs, AcquisitionData)
##        assert isinstance(imgs, ImageData)
        assert_validity(acqs, AcquisitionData)
        assert_validity(imgs, ImageData)
        self.handle = None
        self.handle = \
            pygadgetron.cGT_AcquisitionModel(acqs.handle, imgs.handle)
        check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def set_coil_sensitivity_maps(self, csm):
        '''
        Specifies the coil sensitivity maps to be used by the model.
        csm: CoilSensitivityData
        '''
##        assert isinstance(csm, CoilSensitivityData)
        assert_validity(csm, CoilSensitivityData)
        try_calling(pygadgetron.cGT_setCSMs(self.handle, csm.handle))
    def forward(self, image):
        '''
        Projects an image into (simulated) acquisitions space.
        The resulting acquisition data simulates the actual data
        expected to be received from the scanner.
        image: ImageData
        '''
##        assert isinstance(image, ImageData)
        assert_validity(image, ImageData)
        ad = AcquisitionData()
        ad.handle = pygadgetron.cGT_AcquisitionModelForward\
            (self.handle, image.handle)
        check_status(ad.handle)
        return ad;
    def backward(self, ad):
        '''
        Back-projects acquisition data into image space using a complex
        transpose of the forward projection.
        ad: AcquisitionData
        '''
##        assert isinstance(ad, AcquisitionData)
        assert_validity(ad, AcquisitionData)
        image = ImageData()
        image.handle = pygadgetron.cGT_AcquisitionModelBackward\
            (self.handle, ad.handle)
        check_status(image.handle)
        return image

class Gadget:
    '''
    Class for Gadgetron gadgets.
    '''
    def __init__(self, name):
        '''
        Creates a gadget of specified type and properties.
        name: a string of the form gadget_type(property1=value1, ...)
        '''
        self.handle = None
        name, prop = name_and_parameters(name)
        self.handle = pygadgetron.cGT_newObject(name)
        check_status(self.handle)
        if prop is not None:
            self.set_properties(prop)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def set_property(self, prop, value):
        '''
        Assigns specified value to specified gadget property.
        prop : property name (string)
        value: property value (string)
        '''
        try_calling(pygadgetron.cGT_setGadgetProperty(self.handle, prop, value))
    def set_properties(self, prop):
        '''
        Assigns specified values to specified gadget properties.
        prop: a string with comma-separated list of property value assignments 
              prop_name=prop_value
        '''
        try_calling(pygadgetron.cGT_setGadgetProperties(self.handle, prop))
    def value_of(self, prop):
        '''
        Returns the string representation of the value of specified property.
        prop: property name (string)
        '''
        return _char_par(self.handle, 'gadget', prop)

class GadgetChain:
    '''
    Class for Gadgetron chains.
    '''
    def __init__(self):
        self.handle = pygadgetron.cGT_newObject('GadgetChain')
        check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
##    def add_reader(self, id, reader):
##        '''
##        Adds reader gadget (a gadget that receives data from the client) to the
##        chain.
##        id    : gadget id (string)
##        reader: Gadget of reader type
##        '''
##        assert isinstance(reader, Gadget)
##        try_calling(pygadgetron.cGT_addReader(self.handle, id, reader.handle))
##    def add_writer(self, id, writer):
##        '''
##        Adds writer gadget (a gadget that sends data to the client) to the
##        chain.
##        id    : gadget id (string)
##        writer: Gadget of writer type
##        '''
##        assert isinstance(writer, Gadget)
##        try_calling(pygadgetron.cGT_addWriter(self.handle, id, writer.handle))
    def add_gadget(self, id, gadget):
        '''
        Adds a gadget to the chain.
        id    : gadget id (string)
        writer: Gadget
        '''
        assert isinstance(gadget, Gadget)
        try_calling(pygadgetron.cGT_addGadget(self.handle, id, gadget.handle))
    def set_gadget_property(self, id, prop, value):
        '''
        Assigns specified value to specified gadget property.
        id   : gadget id
        prop : property name (string)
        value: property value (string)
        '''
        if type(value) == type('abc'):
            v = value
        else:
            v = repr(value).lower()
        hg = _parameterHandle(self.handle, 'gadget_chain', id)
        try_calling(pygadgetron.cGT_setGadgetProperty(hg, prop, v))
        pyiutil.deleteDataHandle(hg)
    def value_of_gadget_property(self, id, prop):
        '''
        Returns the string representation of the value of specified property.
        id  : gadget id
        prop: property name (string)
        '''
        hg = _parameterHandle(self.handle, 'gadget_chain', id)
        hv = _parameterHandle(hg, 'gadget', prop)
        value = pyiutil.charDataFromHandle(hv)
        pyiutil.deleteDataHandle(hg)
        pyiutil.deleteDataHandle(hv)
        return value

class Reconstructor(GadgetChain):
    '''
    Class for a chain of gadgets that has AcquisitionData on input and 
    ImageData on output.
    '''
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
        '''
        Sets the input.
        input_data: AcquisitionData
        '''
        assert isinstance(input_data, AcquisitionData)
        self.input_data = input_data
    def process(self):
        '''
        Processes the input with the gadget chain.
        '''
        if self.input_data is None:
            raise error('no input data')
        try_calling(pygadgetron.cGT_reconstructImages\
             (self.handle, self.input_data.handle))
    def get_output(self, subset = None):
        '''
        Returns specified subset of the output ImageData. If no subset is 
        specified, returns all output.
        subset: the name of the subset (e.g. images, gfactors,...)
        '''
        output = ImageData()
        output.handle = pygadgetron.cGT_reconstructedImages(self.handle)
        check_status(output.handle)
        if subset is None:
            return output
        else:
            return output.select('GADGETRON_DataRole', subset)
    def reconstruct(self, input_data):
        '''
        Returns the output from the chain for specified input.
        input_data: AcquisitionData
        '''
##        assert isinstance(input_data, AcquisitionData)
        assert_validity(input_data, AcquisitionData)
        handle = pygadgetron.cGT_reconstructImages\
             (self.handle, input_data.handle)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
        images = ImageData()
        images.handle = pygadgetron.cGT_reconstructedImages(self.handle)
        check_status(images.handle)
        return images

class ImageDataProcessor(GadgetChain):
    '''
    Class for a chain of gadgets that has ImageData on input and output.
    '''
    def __init__(self, list = None):
        '''
        Creates an image processor specified by a list of gadgets.
        list: Python list of gadget description strings, each gadget 
              description being a string of the form
                '[label:]gadget_name[(property1=value1[,...])]'
              (square brackets embrace optional items, ... stands for etc.)
        '''
        self.handle = None
        self.handle = pygadgetron.cGT_newObject('ImagesProcessor')
        check_status(self.handle)
        # TODO: handle input and output in cSTIR
        self.input_data = None
        self.output_data = None
        if list is None:
            return
        for i in range(len(list)):
            label, name = label_and_name(list[i])
            self.add_gadget(label, Gadget(name))
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def set_input(self, input_data):
        '''
        Sets the input data.
        '''
        assert isinstance(input_data, ImageData)
        self.input_data = input_data
    def process(self, input_data = None):
        '''
        Returns the output from the chain.
        input_data: ImageData
        '''
        if input_data is not None:
            self.set_input(input_data)
        if self.input_data is None:
            raise error('input data not set')
##        assert isinstance(self.input_data, ImageData)
        assert_validity(self.input_data, ImageData)
        image = ImageData()
        image.handle = pygadgetron.cGT_processImages\
             (self.handle, self.input_data.handle)
        check_status(image.handle)
        self.output_data = image
        return image
    def get_output(self):
        '''
        Returns the output data.
        '''
        return self.output_data
##    def apply(self, image): # cannot be done this way (Python passes the value)
##        processed_image = self.process(image)
##        image = processed_image

class AcquisitionDataProcessor(GadgetChain):
    '''
    Class for a chain of gadgets that has AcquisitionData on input and output.
    '''
    def __init__(self, list = None):
        '''
        Creates an acquisition processor specified by a list of gadgets.
        list: Python list of gadget description strings, each gadget 
              description being a string of the form
                '[label:]gadget_name[(property1=value1[,...])]'
              (square brackets embrace optional items, ... stands for etc.)
        '''
        self.handle = None
        self.handle = pygadgetron.cGT_newObject('AcquisitionsProcessor')
        check_status(self.handle)
        self.input_data = None
        self.output_data = None
        if list is None:
            return
        for i in range(len(list)):
            label, name = label_and_name(list[i])
            self.add_gadget(label, Gadget(name))
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def set_input(self, input_data):
        '''
        Sets the input data.
        '''
        assert isinstance(input_data, AcquisitionData)
        self.input_data = input_data
    def process(self, input_data = None):
        '''
        Returns the output from the chain for specified input.
        input_data: AcquisitionData
        '''
        if input_data is not None:
            self.set_input(input_data)
        if self.input_data is None:
            raise error('input data not set')
##        assert isinstance(self.input_data, AcquisitionData)
        assert_validity(self.input_data, AcquisitionData)
        acquisitions = AcquisitionData()
        acquisitions.handle = pygadgetron.cGT_processAcquisitions\
             (self.handle, self.input_data.handle)
        check_status(acquisitions.handle)
        self.output_data = acquisitions
        return acquisitions
    def get_output(self):
        '''
        Returns the output data.
        '''
        return self.output_data

class FullySampledReconstructor(Reconstructor):
    '''
    Class for a reconstructor from fully sampled Cartesian raw data.
    '''
    def __init__(self):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject('SimpleReconstructionProcessor')
        check_status(self.handle)
        self.input_data = None
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    
class CartesianGRAPPAReconstructor(Reconstructor):
    '''
    Class for a reconstructor from undersampled Cartesian raw data.
    '''
    def __init__(self):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject\
            ('SimpleGRAPPAReconstructionProcessor')
        check_status(self.handle)
        self.input_data = None
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def compute_gfactors(self, flag):
        self.set_gadget_property('gadget4', 'send_out_gfactor', flag)
    
def preprocess_acquisition_data(input_data):
    '''
    Acquisition processor function that adjusts noise and asymmetrich echo and
    removes readout oversampling.
    '''
    assert isinstance(input_data, AcquisitionData)
    return input_data.process(\
        ['NoiseAdjustGadget', \
         'AsymmetricEchoAdjustROGadget', \
         'RemoveROOversamplingGadget'])
    
