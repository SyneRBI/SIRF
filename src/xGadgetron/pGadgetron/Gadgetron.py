''' 
Object-Oriented wrap for the cGadgetron-to-Python interface pygadgetron.py
'''

# SyneRBI Synergistic Image Reconstruction Framework (SIRF)
# Copyright 2015 - 2020 Rutherford Appleton Laboratory STFC
# Copyright 2018 - 2020 University College London
# Copyright 2018 - 2021 Physikalisch-Technische Bundesanstalt (PTB)
#
# This is software developed for the Collaborative Computational
# Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
# (http://www.ccpsynerbi.ac.uk/).
#
# Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#       http://www.apache.org/licenses/LICENSE-2.0
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

import abc
import numpy
import os
from numbers import Number, Complex, Integral
try:
    import pylab
    HAVE_PYLAB = True
except:
    HAVE_PYLAB = False
import sys
import time
from deprecation import deprecated

from sirf.Utilities import show_2D_array, show_3D_array, error, check_status, \
     try_calling, assert_validity, assert_validities, label_and_name, \
     name_and_parameters, parse_arglist, \
     cpp_int_dtype, \
     examples_data_path, existing_filepath, \
     pTest, RE_PYEXT
import sirf
from sirf import SIRF
from sirf.SIRF import DataContainer
import sirf.pyiutilities as pyiutil
import sirf.pygadgetron as pygadgetron
import sirf.pysirf as pysirf

import sirf.Gadgetron_params as parms

if sys.version_info[0] >= 3 and sys.version_info[1] >= 4:
    ABC = abc.ABC
else:
    ABC = abc.ABCMeta('ABC', (), {})

# max number of acquisitions dimensions
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

# data path finding helper functions
@deprecated(
    deprecated_in="2.0.0", removed_in="4.0", current_version=sirf.__version__,
    details="use examples_data_path() instead")
def mr_data_path():
    '''
    Returns default path to MR raw data files.
    '''
    return examples_data_path('MR')

### low-level client functionality
### likely to be obsolete - not used for a long time
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

class Image(object):
    '''
    Provides access to ISMRMRD::Image parameters (cf. ismrmrd.h).
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
    def version(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'image', 'version')
    def flags(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'image', 'flags')
        #return _int_par(self.handle, 'image', 'flags')
    def data_type(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'image', 'data_type')
    def measurement_uid(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'image', 'measurement_uid')
    def channels(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'image', 'channels')
    def average(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'image', 'average')
    def slice(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'image', 'slice')
    def contrast(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'image', 'contrast')
    def phase(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'image', 'phase')
    def repetition(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'image', 'repetition')
    def set(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'image', 'set')
    def acquisition_time_stamp(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'image', 'acquisition_time_stamp')
    def image_type(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'image', 'image_type')
    def image_index(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'image', 'image_index')
    def image_series_index(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'image', 'image_series_index')
    def attribute_string_len(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'image', 'attribute_string_len')
    def matrix_size(self):
        assert self.handle is not None
        return parms.uint16_pars(self.handle, 'image', 'matrix_size', 3)[::-1]
        #return _uint16_pars(self.handle, 'image', 'matrix_size', 3)[::-1]
    def physiology_time_stamp(self):
        assert self.handle is not None
        return parms.uint32_pars(self.handle, 'image', 'physiology_time_stamp', 3)
    def field_of_view(self):
        assert self.handle is not None
        return parms.float_pars(self.handle, 'image', 'field_of_view', 3)
    def position(self):
        assert self.handle is not None
        return parms.float_pars(self.handle, 'image', 'position', 3)
    def read_dir(self):
        assert self.handle is not None
        return parms.float_pars(self.handle, 'image', 'read_dir', 3)
    def phase_dir(self):
        assert self.handle is not None
        return parms.float_pars(self.handle, 'image', 'phase_dir', 3)
    def slice_dir(self):
        assert self.handle is not None
        return parms.float_pars(self.handle, 'image', 'slice_dir', 3)
    def patient_table_position(self):
        assert self.handle is not None
        #return parms.float_pars \
        return parms.float_pars \
               (self.handle, 'image', 'patient_table_position', 3)
    def info(self, method):
        return eval('self.' + method + '()')

#class ImageData(DataContainer):
class ImageData(SIRF.ImageData):
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

    def abs(self, out=None):
        if out is None:
            images = ImageData()
        else:
            images = out
        images.handle = pygadgetron.cGT_realImageData(self.handle, 'abs')
        check_status(images.handle)
        if out is None:
            return images

    def real(self):
        images = ImageData()
        images.handle = pygadgetron.cGT_realImageData(self.handle, 'real')
        check_status(images.handle)
        return images

    def read_from_file(self, file):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
        self.handle = pygadgetron.cGT_readImages(file)
        check_status(self.handle)

    def from_acquisition_data(self, ad):
        assert isinstance(ad, AcquisitionData), "Please pass a AcquisitionData object"
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
        self.handle = pygadgetron.cGT_ImageFromAcquisitiondata(ad.handle)
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

    def set_ISMRMRD_image_type(self, imtype):
        assert self.handle is not None
        try_calling(pygadgetron.cGT_setImageType(self.handle, imtype))

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

    def image(self, im_num):
        return Image(self, im_num)

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

    def get_ISMRMRD_info(self, par):
        '''
        Returns the array of values of the specified image information 
        parameter. Parameters names are the same as the names of sirf.Gadgetron.Image class
        public methods (except is_real and info).

        par: parameter name (as a string)

        Examples:

        # to get information on the timing:
        slice = image.get_ISMRMRD_info('acquisition_time_stamp')

        # to get the unit vector orthogonal to the slice and directed
        # to the next slice (in ISMRMRD coordinate system):
        slice_dir = image.get_ISMRMRD_info('slice_dir')
        '''
        ni = self.number()
        info = numpy.empty((ni,), dtype = object)
        for i in range(ni):
            image = self.image(i)
            info[i] = image.info(par)
        return info

    @deprecated(details="Please use get_ISMRMRD_info method instead")
    def get_info(self, par):
        return self.get_ISMRMRD_info(par)

    def fill(self, data):
        '''
        Fills self's image data with specified values.
        data: Python Numpy array or ImageData
        '''
        assert self.handle is not None
        
        if isinstance(data, ImageData):
            return super(ImageData, self).fill(data)
        
        if isinstance(data, numpy.ndarray):
            dims = self.dimensions()
            shape = data.shape
            if numpy.prod(shape) != numpy.prod(dims):
                msg = 'cannot fill ImageData of size %s with data of size %s'
                raise ValueError(msg % (repr(dims), repr(shape)))
            the_data = data
            if self.is_real():
                if data.dtype != numpy.float32:
                    the_data = data.astype(numpy.float32)
            else:
                if data.dtype != numpy.complex64:
                    the_data = data.astype(numpy.complex64)
            convert = not data.flags['C_CONTIGUOUS']
            if convert:
                the_data = numpy.ascontiguousarray(the_data)
            if self.is_real():
                try_calling(pygadgetron.cGT_setImageDataFromFloatArray\
                    (self.handle, the_data.ctypes.data))
            else:
                try_calling(pygadgetron.cGT_setImageDataFromCmplxArray\
                    (self.handle, the_data.ctypes.data))
        elif isinstance (data, Complex):
            arr = data + numpy.zeros(self.shape, dtype=numpy.complex64)
            return self.fill(arr)
        elif isinstance (data, Number):
            arr = data + numpy.zeros(self.shape, dtype=numpy.float32)
            return self.fill(arr)
        else:
            raise error('wrong fill value.' + \
                        ' Should be ImageData, numpy.ndarray or number. Got {}'.format(type(data)))
        return self

    def dimensions(self):
        '''
        Returns the dimensions of 3D/4D Numpy ndarray of all self's images.
        '''
        if self.number() < 1:
            return 0
        assert self.handle is not None
        dim = numpy.ndarray((4,), dtype=cpp_int_dtype())
        image = Image(self)
        pygadgetron.cGT_getImageDim(image.handle, dim.ctypes.data)
        nx = dim[0]
        ny = dim[1]
        nz = dim[2]
        nc = dim[3]
        
        nz = nz*self.number()
        
        if nc == 1: # for backward compatibility
            return nz, ny, nx
        else:
            return nc, nz, ny, nx
        
    def as_array(self):
        '''
        Returns all self's images as a 3D or 4D Numpy ndarray.
        '''
        dims = self.dimensions()
        
        assert self.handle is not None
        if self.number() < 1:
            return numpy.ndarray((0,), dtype = numpy.float32)
        if self.is_real():
            array = numpy.ndarray(dims, dtype = numpy.float32)
            try_calling(pygadgetron.cGT_getImageDataAsFloatArray\
                (self.handle, array.ctypes.data))
        else:
            array = numpy.ndarray(dims, dtype = numpy.complex64)
            try_calling(pygadgetron.cGT_getImageDataAsCmplxArray\
                (self.handle, array.ctypes.data))
                
        if len(dims) != 4:
            return array

        nc, nz, ny, nx = dims
        ns = self.number() # number of total dynamics (slices, contrasts, etc.)
        nz = nz//ns        # z-dimension of a slice

        # hope Numpy is clever enough to do all this in-place:
        array = numpy.reshape(array, (ns, nc, nz, ny, nx))
        array = numpy.swapaxes(array, 0, 1)
        array = numpy.reshape(array, (nc, ns*nz, ny, nx))

        return array
                
    def copy(self):
        '''alias of clone'''
        return self.clone()
    def show(self, zyx=None, slice=None, title=None, cmap='gray', postpone=False):
        '''Displays xy-cross-section(s) of images.'''
        assert self.handle is not None
        if not HAVE_PYLAB:
            print('pylab not found')
            return
        data = self.as_array()
        nz = data.shape[0]
        if isinstance(slice, (Integral,numpy.integer)):
            if slice < 0 or slice >= nz:
                return
            ni = 1
            slice = [slice]
        elif slice is None:
            ni = nz
            slice = range(nz)
        else:
            try:
                ni = len(slice)
            except:
                raise error('wrong slice list')
        if title is None:
            title = 'Selected images'
        if ni >= 16:
            tiles = (4, 4)
        else:
            tiles = None
        f = 0
        while f < ni:
            t = min(f + 16, ni)
            err = show_3D_array(abs(data), index=slice[f : t], \
                                tile_shape=tiles, cmap=cmap, \
                                zyx=zyx, label='image', \
                                xlabel='samples', ylabel='readouts', \
                                suptitle=title, \
                                show=(t == ni) and not postpone)
            f = t
    def allocate(self, value=0, **kwargs):
        '''Method to allocate an ImageData and set its values
        
        CIL/SIRF compatibility
        '''
        if value in ['random', 'random_int']:
            out = self.clone()
            shape = out.as_array().shape
            seed = kwargs.get('seed', None)
            if seed is not None:
                numpy.random.seed(seed) 
            if value == 'random':
                out.fill(numpy.random.random_sample(shape))
            elif value == 'random_int':
                max_value = kwargs.get('max_value', 100)
                out.fill(numpy.random.randint(max_value,size=shape))
        else:
            out = self.clone()
            tmp = value * numpy.ones(out.as_array().shape)
            out.fill(tmp)
        return out

    def print_header(self, im_num):
        """Print the header of one of the images. zero based."""
        try_calling(pygadgetron.cGT_print_header(self.handle, im_num))
    @property
    def dtype(self):
        if self.is_real():
            return numpy.float32
        return numpy.complex64
    @property
    def shape(self):
        return self.dimensions()
        
SIRF.ImageData.register(ImageData)


class CoilImagesData(ImageData):
    '''
    Class for a coil images (ci) container.
    Each item in the container is a 4D complex array of coil images values
    on an xyz-slice.
    '''
    def __init__(self):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject('CoilImages')
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def same_object(self):
        return CoilImagesData()
    def calculate(self, acq):
        dcw = compute_kspace_density(acq)
        acq = acq * dcw
        try_calling(pygadgetron.cGT_computeCoilImages(self.handle, acq.handle))

SIRF.ImageData.register(CoilImagesData)


class CoilSensitivityData(ImageData):
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
    def calculate(self, data, method=None):
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
        nit = self.smoothness
        
        if method is not None:
            method_name, parm_list = name_and_parameters(method)
            parm = parse_arglist(parm_list)
            if 'niter' in parm:
                nit = int(parm['niter'])
        else:
            method_name = 'SRSS'
            parm = {}
        
        parms.set_int_par(self.handle, 'coil_sensitivity', 'smoothness', nit)

        if isinstance(data, AcquisitionData):
            self.__calc_from_acquisitions(data, method_name)
        elif isinstance(data, CoilImagesData):
            self.__calc_from_images(data, method_name)
        else:
            raise error('Cannot calculate coil sensitivities from %s' % \
                        repr(type(data)))

    def __calc_from_acquisitions(self, data, method_name):

        if data.handle is None:
            raise AssertionError("The handle for data is None. Please pass valid acquisition data.")

        dcw = compute_kspace_density(data)

        data = data * dcw
        if method_name == 'Inati':
            try:
                from ismrmrdtools import coils
            except:
                raise error('Inati method requires ismrmrd-python-tools')

            cis = CoilImagesData()
            try_calling(pygadgetron.cGT_computeCoilImages(cis.handle, data.handle))
            cis_array = cis.as_array()
            csm, _ = coils.calculate_csm_inati_iter(cis_array)

            if self.handle is not None:
                pyiutil.deleteDataHandle(self.handle)
            self.handle = pysirf.cSIRF_clone(cis.handle)
            nc, nz, ny, nx = self.dimensions()
            ns = self.number() # number of total dynamics (slices, contrasts, etc.)
            nz = nz//ns        # z-dimension of a slice
            csm = numpy.reshape(csm, (nc, ns, nz, ny, nx))
            csm = numpy.swapaxes(csm, 0,  1)
            self.fill(csm.astype(numpy.complex64))
        
        elif method_name == 'SRSS':
            try_calling(pygadgetron.cGT_computeCoilSensitivities(self.handle, data.handle))

    def __calc_from_images(self, data, method_name):

        if data.handle is None:
            raise AssertionError("The handle for data is None. Please pass valid image data.")


        if method_name == 'Inati':
            try:
                from ismrmrdtools import coils
            except:
                raise error('Inati method requires ismrmrd-python-tools')
                
            cis_array = data.as_array()
            csm, _ = coils.calculate_csm_inati_iter(cis_array)
            if self.handle is not None:
                pyiutil.deleteDataHandle(self.handle)
            self.handle = pysirf.cSIRF_clone(data.handle)

            nc, nz, ny, nx = self.dimensions()
            ns = self.number() # number of total dynamics (slices, contrasts, etc.)
            nz = nz//ns        # z-dimension of a slice
            csm = numpy.reshape(csm, (nc, ns, nz, ny, nx))
            csm = numpy.swapaxes(csm, 0,  1)

            self.fill(csm.astype(numpy.complex64))

        elif method_name == 'SRSS':

            try_calling(pygadgetron.cGT_computeCoilSensitivitiesFromCoilImages \
                (self.handle, data.handle))

        else:
            raise error('Unknown method %s' % method_name)

    def __calc_from_acquisitions(self, data, method_name):
        assert data.handle is not None
        dcw = compute_kspace_density(data)
        data = data * dcw
        cis = CoilImagesData()
        try_calling(pygadgetron.cGT_computeCoilImages(cis.handle, data.handle))
        self.__calc_from_images(cis, method_name)

DataContainer.register(CoilSensitivityData)


class Acquisition(object):
    ''' Provides access to ISMRMRD::Acquisition parameters (cf. ismrmrd.h).
    '''
    def __init__(self, file = None):
        self.handle = None
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def version(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'acquisition', 'version')
    def flags(self):
        '''
        Returns acquisition flags as an integer (each bit corresponding to a 
        flag).
        '''
        assert self.handle is not None
        return parms.int_par(self.handle, 'acquisition', 'flags')
        #return _int_par(self.handle, 'acquisition', 'flags')
    def measurement_uid(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'acquisition', 'measurement_uid')
    def scan_counter(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'acquisition', 'scan_counter')
    def acquisition_time_stamp(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'acquisition', 'acquisition_time_stamp')
    def number_of_samples(self):
        '''
        returns the number of samples in the readout direction.
        '''
        assert self.handle is not None
        return parms.int_par(self.handle, 'acquisition', 'number_of_samples')
    def available_channels(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'acquisition', 'available_channels')
    def active_channels(self):
        '''
        Returns the number of active channels (coils).
        '''
        assert self.handle is not None
        return parms.int_par(self.handle, 'acquisition', 'active_channels')
    def discard_pre(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'acquisition', 'discard_pre')
    def discard_post(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'acquisition', 'discard_post')
    def center_sample(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'acquisition', 'center_sample')
    def encoding_space_ref(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'acquisition', 'encoding_space_ref')
    def trajectory_dimensions(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'acquisition', 'trajectory_dimensions')
    def kspace_encode_step_1(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'acquisition', 'idx_kspace_encode_step_1')
    def kspace_encode_step_2(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'acquisition', 'idx_kspace_encode_step_2')
    def average(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'acquisition', 'idx_average')
    def slice(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'acquisition', 'idx_slice')
    def contrast(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'acquisition', 'idx_contrast')
    def phase(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'acquisition', 'idx_phase')
    def repetition(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'acquisition', 'idx_repetition')
    def set(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'acquisition', 'idx_set')
    def segment(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'acquisition', 'idx_segment')
    def physiology_time_stamp(self):
        assert self.handle is not None
        return parms.uint32_pars(self.handle, 'acquisition', 'physiology_time_stamp', 3)
    def channel_mask(self):
        assert self.handle is not None
        return parms.uint64_pars(self.handle, 'acquisition', 'channel_mask', 16)
    def sample_time_us(self):
        assert self.handle is not None
        return parms.float_par(self.handle, 'acquisition', 'sample_time_us')
    def position(self):
        assert self.handle is not None
        return parms.float_pars(self.handle, 'acquisition', 'position', 3)
    def read_dir(self):
        assert self.handle is not None
        return parms.float_pars(self.handle, 'acquisition', 'read_dir', 3)
    def phase_dir(self):
        assert self.handle is not None
        return parms.float_pars(self.handle, 'acquisition', 'phase_dir', 3)
    def slice_dir(self):
        assert self.handle is not None
        return parms.float_pars(self.handle, 'acquisition', 'slice_dir', 3)
    def patient_table_position(self):
        assert self.handle is not None
        return parms.float_pars \
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
        '''Sets acquisition data storage scheme.'''
        if scheme != 'memory':
            msg = 'WARNING: storage scheme ' + repr(scheme) + ' not supported,'
            msg += ' using memory storage scheme instead'
            print(msg)
    @staticmethod
    def get_storage_scheme():
        '''Returns acquisition data storage scheme.
        '''
        return 'memory'
    def same_object(self):
        return AcquisitionData()
    def new_acquisition_data(self, empty=True):
        new_ad = AcquisitionData()
        if empty:
            new_ad.handle = pygadgetron.cGT_createEmptyAcquisitionData(self.handle)
        else:
            new_ad.handle = pygadgetron.cGT_cloneAcquisitions(self.handle)
        check_status(new_ad.handle)
        return new_ad
    def number_of_readouts(self, select='image'):
        if select == 'image':
            dim = self.dimensions()
            return dim[0]
        else:
            return self.number()
    def number_of_acquisitions(self, select='image'):
        return self.number_of_readouts

    def check_traj_type(self, trajname):
        '''
        Checks if the data is of the trajectory type trajname.
        trajname: string with trajectory name.
        Possible choices are:
            - cartesian
            - radial
            - epi
            - goldenangle
            - spiral
            - other
        '''
        list_available_trajs = ('cartesian', 'epi', 'radial', 'goldenangle', 'spiral', 'other')
        if trajname not in list_available_trajs:
            raise AssertionError("The trajectory you asked for is not among the available trajectoryies")

        xml_hdr = self.get_header()
        traj_id_substring = "<trajectory>"+trajname+"</trajectory>"

        return traj_id_substring in xml_hdr
     
    def sort(self):
        '''
        Sorts acquisitions (currently, with respect to acquisition_time_stamp)
        '''
        assert self.handle is not None
        try_calling(pygadgetron.cGT_sortAcquisitions(self.handle))
        self.sorted = True
    def sort_by_time(self):
        '''
        Sorts acquisitions with respect to acquisition_time_stamp
        '''
        assert self.handle is not None
        try_calling(pygadgetron.cGT_sortAcquisitionsByTime(self.handle))
    def is_sorted(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'acquisitions', 'sorted')
        #return self.sorted
    def is_undersampled(self):
        assert self.handle is not None
        return parms.int_par(self.handle, 'acquisitions', 'undersampled')
    def set_header(self, header):
        assert self.handle is not None
        try_calling(pygadgetron.cGT_setAcquisitionsInfo(self.handle, header))
    def get_header(self):
        assert self.handle is not None
        return parms.char_par(self.handle, 'acquisitions', 'info')
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
    def acquisition(self, num):
        '''
        Returns the specified acquisition.
        num: acquisition number
        '''
        assert self.handle is not None
        acq = Acquisition()
        acq.handle = pygadgetron.cGT_acquisitionFromContainer(self.handle, int(num))
        check_status(acq.handle)
        return acq
    def append_acquisition(self, acq):
        '''
        Appends acquistion to AcquisitionData.
        '''
        assert self.handle is not None
        try_calling( pygadgetron.cGT_appendAcquisition(self.handle, acq.handle))
    
    def get_subset(self, idx):
        '''
        Returns AcquisitionData object with subset of acquisitions defined by idx
        '''
        assert self.handle is not None
        subset = AcquisitionData()
        idx = numpy.array(idx, dtype = cpp_int_dtype())
        subset.handle = pygadgetron.cGT_getAcquisitionsSubset(self.handle, idx.ctypes.data, idx.size)
        check_status(subset.handle)
        
        return subset
    
    def set_user_floats(self, data, idx):
        '''
        Writes the data into the user_float[idx] data field of the acquisition
        data header of each acquisition in the container to pass additional data
        into the raw data.
        data: numpy array
        idx: integer in range 0 to 7
        '''        
        if self.handle is None:
            raise AssertionError('self.handle is None')
                    
        if data.size != self.number():
            raise AssertionError('Please give as many datapoints as there are acquisitions')
        
        if idx > 7 or idx < 0 or not isinstance(idx,int):
            raise AssertionError('Please give an integer from [0,...,7]')

        if data.dtype is not numpy.float32:
            data = data.astype(numpy.float32)

        try_calling(pygadgetron.cGT_setAcquisitionUserFloat\
                    (self.handle, data.ctypes.data, idx))

    def dimensions(self):
        '''
        Returns acquisitions dimensions as a tuple (na, nc, ns), where na is
        the number of acquisitions, nc the number of coils and ns the number of
        samples.
        '''
        assert self.handle is not None
        if self.number() < 1:
            return numpy.zeros((MAX_ACQ_DIMENSIONS,), dtype=cpp_int_dtype())
        dim = numpy.ones((MAX_ACQ_DIMENSIONS,), dtype=cpp_int_dtype())
        hv = pygadgetron.cGT_getAcquisitionDataDimensions\
             (self.handle, dim.ctypes.data)
        pyiutil.deleteDataHandle(hv)
        dim[2] = numpy.prod(dim[2:])
        return tuple(dim[2::-1])

    def get_ISMRMRD_info(self, par, which='all'):
        '''
        Returns the array of values of the specified acquisition information 
        parameter.

        par: parameter name (see sirf.Gadgetron.Acquisition class methods except info)
        which: specifies the range of acquisitions whose parameters are returned

        Example:
        # to retrieve readouts flags for acquisitions 0 to 10:
        flags = acq_data.get_ISMRMRD_info('flags', range(10))

        # for phase encoding information
        encoding = acq_data.get_ISMRMRD_info('kspace_encode_step_1')

        '''
        na, nc, ns = self.dimensions()
        if which == 'all':
            rng = range(na)
        else:
            rng = which
            na = len(rng)
        f = min(rng)
        t = max(rng) + 1
        info = numpy.ndarray((2,), dtype=cpp_int_dtype())
        try_calling(pygadgetron.cGT_acquisitionParameterInfo \
                    (self.handle, par, info.ctypes.data))
        n = int(info[1])
        if info[0] == 0:
            values = numpy.ndarray((na, n), dtype=numpy.uint64)
            try_calling(pygadgetron.cGT_acquisitionParameterValuesInt \
                        (self.handle, par, f, t, n, values.ctypes.data))
        else:
            values = numpy.ndarray((na, n), dtype=numpy.float32)
            try_calling(pygadgetron.cGT_acquisitionParameterValuesFloat \
                        (self.handle, par, f, t, n, values.ctypes.data))
        if n == 1:
            values = numpy.reshape(values, (na,))
        return values

        # Python way is much slower for large data
        info = numpy.empty((na,), dtype = object)
        i = 0
        for a in rng:
            acq = self.acquisition(a)
            info[i] = acq.info(par)
            i += 1
        return info

    @deprecated(details="Please use the get_ISMRMRD_info method instead")
    def get_info(self, par, which='all'):
        return self.get_ISMRMRD_info(par, which)

    def fill(self, data, select='image'):
        '''
        Fills self's acquisitions with specified values.
        data: Python Numpy array or AcquisitionData
        select: specifies whether all or only image-related acquisitions are
                filled with values from the array data.
        '''
        assert self.handle is not None
        if isinstance(data, AcquisitionData):
            try_calling(pygadgetron.cGT_fillAcquisitionDataFromAcquisitionData\
                (self.handle, data.handle))
            return
        elif isinstance(data, numpy.ndarray):
            dims = self.dimensions()
            shape = data.shape
            if shape != dims:
                msg = 'cannot fill AcquisitionData of size %s' \
                      + ' with data of size %s'
                raise ValueError(msg % (repr(dims), repr(shape)))
            if data.dtype is not numpy.complex64:
                the_data = data.astype(numpy.complex64)
            else:
                the_data = data
            convert = not data.flags['C_CONTIGUOUS']
            if convert:
                the_data = numpy.ascontiguousarray(the_data)
            if select == 'all':
                fill_all = 1
            else: # fill only image-related
                fill_all = 0
            try_calling(pygadgetron.cGT_fillAcquisitionData\
                (self.handle, the_data.ctypes.data, fill_all))
        elif isinstance (data, Complex):
            arr = data + numpy.zeros(self.shape, dtype=numpy.complex64)
            return self.fill(arr)
        elif isinstance (data, Number):
            arr = data + numpy.zeros(self.shape, dtype=numpy.float32)
            return self.fill(arr)
        else:
            raise error('wrong fill value.' + \
                        ' Should be AcquisitionData, numpy.ndarray or number. Got {}'.format(type(data)))
        return self

    def as_array(self, acq=None):
        '''
        Returns selected self's acquisition(s) data as a 2D or 3D Numpy ndarray.
        acq: acquisition number; if None all acquisitions data is returned.
        '''
        assert self.handle is not None
        if acq is None:
            ''' return 3D array of all acquisition data
            '''
            na, nc, ns = self.dimensions()
            z = numpy.ndarray((na, nc, ns), dtype = numpy.complex64)
            acq = -1
        else:
            ''' return 2D array of the specified acquisition data
            '''
            a = self.acquisition(acq)
            nc = a.active_channels()
            ns = a.number_of_samples()
            z = numpy.ndarray((nc, ns), dtype = numpy.complex64)            
        try_calling(pygadgetron.cGT_acquisitionDataAsArray\
            (self.handle, z.ctypes.data, acq))
        return z

    def show(self, slice=None, title=None, cmap='gray', power=0.2, \
             postpone=False):
        '''Displays xy-cross-section(s) of images.'''
        assert self.handle is not None
        if not HAVE_PYLAB:
            print('pylab not found')
            return
        data = numpy.transpose(self.as_array(), (1, 0, 2))
        nz = data.shape[0]
        if isinstance(slice, (Integral,numpy.integer)):
            if slice < 0 or slice >= nz:
                return
            ns = 1
            slice = [slice]
##            show_2D_array('slice %d' % slice, data[slice,:,:])
##            return
        elif slice is None:
            ns = nz
            slice = range(nz)
        else:
            try:
                ns = len(slice)
            except:
                raise error('wrong slice list')
        if title is None:
            title = 'Selected images'
        if ns >= 16:
            tiles = (4, 4)
        else:
            tiles = None
        f = 0
        while f < ns:
            t = min(f + 16, ns)
            err = show_3D_array(abs(data), index = slice[f : t], \
                                tile_shape = tiles, \
                                label = 'coil', xlabel = 'samples', \
                                ylabel = 'readouts', \
                                suptitle = title, cmap = cmap, power = power, \
                                show = (t == ns) and not postpone)
            f = t
    
    def allocate(self, value=0, **kwargs):
        '''Method to allocate an AcquisitionData and set its values
        
        CIL/SIRF compatibility
        '''
        if value in ['random', 'random_int']:
            out = self.clone()
            shape = out.as_array().shape
            seed = kwargs.get('seed', None)
            if seed is not None:
                numpy.random.seed(seed) 
            if value == 'random':
                out.fill(numpy.random.random_sample(shape))
            elif value == 'random_int':
                max_value = kwargs.get('max_value', 100)
                out.fill(numpy.random.randint(max_value,size=shape))
        else:
            out = self.clone()
            tmp = value * numpy.ones(out.as_array().shape)
            out.fill(tmp)
        return out
    @property
    def shape(self):
        return self.dimensions()
    @property
    def dtype(self):
        return numpy.complex64
    
    
DataContainer.register(AcquisitionData)

class AcquisitionModel(object):
    '''
    Class for MR acquisition model, an operator that maps images into
    simulated acquisitions.
    '''
    def __init__(self, acqs = None, imgs = None):
        self.handle = None
        if acqs == None:
            self.handle = pygadgetron.cGT_newObject('AcquisitionModel')
        else:
            assert_validity(acqs, AcquisitionData)
            assert_validity(imgs, ImageData)
            self.handle = \
                pygadgetron.cGT_AcquisitionModel(acqs.handle, imgs.handle)
        check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def set_up(self, acqs, imgs):
        assert_validity(acqs, AcquisitionData)
        assert_validity(imgs, ImageData)
        try_calling(pygadgetron.cGT_setUpAcquisitionModel \
            (self.handle, acqs.handle, imgs.handle))
    def set_coil_sensitivity_maps(self, csm):
        '''
        Specifies the coil sensitivity maps to be used by the model.
        csm: CoilSensitivityData
        '''
        assert_validity(csm, CoilSensitivityData)
        try_calling(pygadgetron.cGT_setAcquisitionModelParameter \
            (self.handle, 'coil_sensitivity_maps', csm.handle))
    def norm(self, num_iter=2, verb=0):
        '''Computes the norm of the forward projection operator.
        '''
        assert self.handle is not None
        handle = pygadgetron.cGT_acquisitionModelNorm(self.handle, num_iter, verb)
        check_status(handle)
        r = pyiutil.floatDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return r;
    def forward(self, image):
        '''
        Projects an image into (simulated) acquisitions space.
        The resulting acquisition data simulates the actual data
        expected to be received from the scanner.
        image: ImageData
        '''
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
        assert_validity(ad, AcquisitionData)
        image = ImageData()
        image.handle = pygadgetron.cGT_AcquisitionModelBackward\
            (self.handle, ad.handle)
        check_status(image.handle)
        return image
    def inverse(self, ad, dcw=None):
        '''
        Weights acquisition data with k-space density prior to back-projection
        into image space using a complex transpose of the forward projection.
        ad: AcquisitionData
        dcw: AcquisitionData
        '''
        assert_validity(ad, AcquisitionData)
        
        if dcw is not None:
            assert_validity(dcw, AcquisitionData)
            if ad.shape != dcw.shape:
                raise AssertionError("The shape of the density weights and the acquisition data must be the same.")

        if dcw is None:
            dcw = compute_kspace_density(ad)
            
        ad = ad * dcw

        image = ImageData()
        image.handle = pygadgetron.cGT_AcquisitionModelBackward\
            (self.handle, ad.handle)
        check_status(image.handle)
        return image

    def direct(self, image, out = None):
        '''Alias of forward

           Added for CCPi CIL compatibility
           https://github.com/CCPPETMR/SIRF/pull/237#issuecomment-439894266
        '''
        if out is not None:
            #raise error('out is not supported')
            tmp = self.forward(image)
            out.fill(tmp)
            return
        return self.forward(image)
    def adjoint(self, ad , out = None):
        '''Alias of backward

           Added for CCPi CIL compatibility
           https://github.com/CCPPETMR/SIRF/pull/237#issuecomment-439894266
        '''
        if out is not None:
            #raise error('out is not supported')
            tmp = self.backward(ad)
            out.fill(tmp)
            return
        return self.backward(ad)
    def is_affine(self):
        '''Returns if the acquisition model is affine (i.e. corresponding to A*x+b)'''
        return True
    def is_linear(self):
        '''Returns whether the acquisition model is linear (i.e. corresponding to A*x, with zero background term)'''
        return True

    def range_geometry(self):
        '''Returns the template of AcquisitionData'''
        geom = AcquisitionData()
        geom.handle = pygadgetron.cGT_parameter(
            self.handle, 'AcquisitionModel', 'range geometry')
        check_status(geom.handle)
        return geom

    def domain_geometry(self):
        '''Returns the template of ImageData'''
        geom = ImageData()
        geom.handle = pygadgetron.cGT_parameter(
            self.handle, 'AcquisitionModel', 'domain geometry')
        check_status(geom.handle)
        return geom

class Gadget(object):
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
        return parms.char_par(self.handle, 'gadget', prop)

class GadgetChain(object):
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
    def set_host(self, host):
        '''
        Sets Gadgetron server host.
        host : host name (string)
        '''
        try_calling(pygadgetron.cGT_setHost(self.handle, host))
    def set_port(self, port):
        '''
        Sets Gadgetron server port.
        port : port number (as a string)
        '''
        try_calling(pygadgetron.cGT_setPort(self.handle, port))
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
        hg = parms.parameter_handle(self.handle, 'gadget_chain', id)
        try_calling(pygadgetron.cGT_setGadgetProperty(hg, prop, v))
        pyiutil.deleteDataHandle(hg)
    def value_of_gadget_property(self, id, prop):
        '''
        Returns the string representation of the value of specified property.
        id  : gadget id
        prop: property name (string)
        '''
        hg = parms.parameter_handle(self.handle, 'gadget_chain', id)
        hv = parms.parameter_handle(hg, 'gadget', prop)
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
        self.handle = pygadgetron.cGT_newObject('SimpleReconstructionprocessor')
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
            ('SimpleGRAPPAReconstructionprocessor')
        check_status(self.handle)
        self.input_data = None
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def compute_gfactors(self, flag):
        self.set_gadget_property('gadget4', 'send_out_gfactor', flag)
    
def preprocess_acquisition_data(input_data):
    '''
    Acquisition processor function that adjusts noise and asymmetric echo and
    removes readout oversampling.
    '''
    assert_validity(input_data, AcquisitionData)
    return input_data.process(\
        ['NoiseAdjustGadget', \
         'AsymmetricEchoAdjustROGadget', \
         'RemoveROOversamplingGadget'])
    
def set_grpe_trajectory(ad):
    '''
    Function that fills the trajectory of AcquisitionData with golden angle radial
    phase encoding trajectory.
    ad: AcquisitionData
    '''    
    assert_validity(ad, AcquisitionData)

    try_calling(pygadgetron.cGT_setGRPETrajectory(ad.handle))
    return ad
    
def set_radial2D_trajectory(ad):
    '''
    Function that fills the trajectory of AcquisitionData with linear increment 2D radial
    readout trajectory.
    ad: AcquisitionData
    '''
    assert_validity(ad, AcquisitionData)

    try_calling(pygadgetron.cGT_setRadial2DTrajectory(ad.handle))
    return ad

def set_goldenangle2D_trajectory(ad):
    '''
    Function that fills the trajectory of AcquisitionData with golden angle increment 2D radial
    readout trajectory.
    ad: AcquisitionData
    '''
    assert_validity(ad, AcquisitionData)

    try_calling(pygadgetron.cGT_setGoldenAngle2DTrajectory(ad.handle))
    return ad

def get_data_trajectory(ad):
    '''
    Function that gets the trajectory of AcquisitionData depending on the rawdata trajectory.
    ad: AcquisitionData
    '''    
    assert_validity(ad, AcquisitionData)
    
    if ad.check_traj_type('cartesian'):
        num_traj_pts = ad.number()
        traj_dim = 2
    elif ad.check_traj_type('other'):
        num_traj_pts = ad.number()
        traj_dim = 3
    elif ad.check_traj_type('radial') or ad.check_traj_type('goldenangle'):
        num_traj_pts = ad.number() * ad.dimensions()[2]
        traj_dim = 2
        
    dims = (num_traj_pts, traj_dim)
    traj = numpy.ndarray(dims, dtype = numpy.float32)
    
    try_calling(pygadgetron.cGT_getDataTrajectory(ad.handle, traj.ctypes.data))
    
    return traj

def compute_kspace_density(ad):
    '''
    Function that computes the kspace density depending the
    ad: AcquisitionData
    '''  
    assert_validity(ad, AcquisitionData)

    if ad.check_traj_type('cartesian'):
        return calc_cartesian_dcw(ad)
    elif ad.check_traj_type('other'):
        return calc_rpe_dcw(ad)
    elif ad.check_traj_type('radial') or ad.check_traj_type('goldenangle'):
        return calc_radial_dcw(ad)
    	
    else:
        raise AssertionError("Please only try to recon trajectory types cartesian or other")
    

def calc_cartesian_dcw(ad):
    '''
    Function that computes the kspace weight for a cartesian acquisition by
    averaging out phase encoding points that were acquired multiple times.
    ad: AcquisitionData
    '''
    traj = numpy.transpose(get_data_trajectory(ad))
    traj, inverse, counts = numpy.unique(traj, return_inverse=True, return_counts=True, axis=1)
    
    density_weight = (1.0 / counts)[inverse]
    
    density_weight = numpy.expand_dims(density_weight, axis=1)
    density_weight = numpy.expand_dims(density_weight, axis=2)
    density_weight = numpy.tile(density_weight, (1, ad.shape[1], ad.shape[2]))
    
    dcw = ad.copy()
    dcw.fill(density_weight)
    
    return dcw

def calc_rpe_dcw(ad):
    '''
    Function that computes the kspace weight depending on the distance to the center
    as in a filtered back-projection. Strictly valid only for equally angular-spaced
    radially distributed points
    ad: AcquisitionData
    '''

    traj = numpy.transpose(get_data_trajectory(ad))
    traj = traj[1:3,:]
    ramp_filter = numpy.linalg.norm(traj, axis=0)

    traj, inverse, counts = numpy.unique(traj, return_inverse=True, return_counts=True, axis=1)
    
    num_angles = numpy.max(counts)
    
    density_weight = ( 1.0 / counts)[inverse]  + num_angles * ramp_filter
    
    max_traj_rad = numpy.max(numpy.linalg.norm(traj, axis=0))
    density_weight_norm =  numpy.sum(density_weight) / (max_traj_rad**2 * numpy.pi)
    density_weight = density_weight / density_weight_norm

    density_weight = numpy.expand_dims(density_weight, axis=(1,2))
    density_weight = numpy.tile(density_weight, (1, ad.shape[1], ad.shape[2]))
    
    dcw = ad.copy()
    dcw.fill(density_weight)
    
    return dcw



def calc_radial_dcw(ad):
    '''
    Function that computes the kspace weight depending on the distance to the center
    as in a filtered back-projection. Stricly valid only for equally angular-spaced
    radially distributed points
    ad: AcquisitionData
    '''

    traj = numpy.transpose(get_data_trajectory(ad))
    (na, nc, ns) = ad.dimensions()
  
    ramp_filter = numpy.linalg.norm(traj, axis=0)
    traj, inverse, counts = numpy.unique(traj, return_inverse=True, return_counts=True, axis=1)
    num_angles = numpy.max(counts)
    
    density_weight = ( 1.0 / counts)[inverse]  + num_angles * ramp_filter
    
    max_traj_rad = numpy.max(numpy.linalg.norm(traj, axis=0))
    density_weight_norm =  numpy.sum(density_weight) / (max_traj_rad**2 * numpy.pi)
    density_weight = density_weight / density_weight_norm

    density_weight = numpy.transpose(density_weight)
    density_weight = numpy.expand_dims(density_weight, axis=(1,2))
    density_weight = numpy.reshape(density_weight, (na, 1, ns))
    density_weight = numpy.tile(density_weight, (1, nc, 1))
    
    dcw = ad.copy()
    dcw.fill(density_weight)
    
    return dcw
