''' 
Object-Oriented wrap for the cSTIR-to-Python interface pystir.py
'''

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
## Copyright 2015 - 2017 University College London
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
import inspect
import numpy
import os
try:
    import pylab
    HAVE_PYLAB = True
except:
    HAVE_PYLAB = False
import sys
import time

from pUtilities import *
import pyiutilities as pyiutil
import pystir

try:
    input = raw_input
except NameError:
    pass

if sys.version_info[0] >= 3 and sys.version_info[1] >= 4:
    ABC = abc.ABC
else:
    ABC = abc.ABCMeta('ABC', (), {})

INFO_CHANNEL = 0
WARNING_CHANNEL = 1
ERROR_CHANNEL = 2
ALL_CHANNELS = -1

###########################################################
############ Utilities for internal use only ##############
def _setParameter(hs, set, par, hv, stack = None):
##    try_calling(pystir.cSTIR_setParameter(hs, set, par, hv))
    if stack is None:
        stack = inspect.stack()[1]
    h = pystir.cSTIR_setParameter(hs, set, par, hv)
    check_status(h, stack)
    pyiutil.deleteDataHandle(h)
def _set_char_par(handle, set, par, value):
    h = pyiutil.charDataHandle(value)
    _setParameter(handle, set, par, h, inspect.stack()[1])
    pyiutil.deleteDataHandle(h)
def _set_int_par(handle, set, par, value):
    h = pyiutil.intDataHandle(value)
    _setParameter(handle, set, par, h, inspect.stack()[1])
    pyiutil.deleteDataHandle(h)
def _set_float_par(handle, set, par, value):
    h = pyiutil.floatDataHandle(value)
    _setParameter(handle, set, par, h, inspect.stack()[1])
    pyiutil.deleteDataHandle(h)
def _char_par(handle, set, par):
    h = pystir.cSTIR_parameter(handle, set, par)
    check_status(h, inspect.stack()[1])
    value = pyiutil.charDataFromHandle(h)
    pyiutil.deleteDataHandle(h)
    return value
def _int_par(handle, set, par):
    h = pystir.cSTIR_parameter(handle, set, par)
    check_status(h, inspect.stack()[1])
    value = pyiutil.intDataFromHandle(h)
    pyiutil.deleteDataHandle(h)
    return value
def _float_par(handle, set, par):
    h = pystir.cSTIR_parameter(handle, set, par)
    check_status(h, inspect.stack()[1])
    value = pyiutil.floatDataFromHandle(h)
    pyiutil.deleteDataHandle(h)
    return value
def _getParameterHandle(hs, set, par):
    handle = pystir.cSTIR_parameter(hs, set, par)
    check_status(handle, inspect.stack()[1])
    return handle
def _tmp_filename():
    return repr(int(1000*time.time()))
###########################################################

class MessageRedirector:
    '''
    Class for STIR printing redirection to files/stdout/stderr.
    '''
    def __init__(self, info = None, warn = 'stdout', errr = 'stdout'):
        '''
        Creates MessageRedirector object that redirects STIR's ouput
        produced by info(), warning() and error(0 functions to destinations
        specified respectively by info, warn and err arguments.
        The argument values other than None, stdout, stderr, cout and cerr
        are interpreted as filenames.
        None and empty string value suppresses printing.
        '''
        if info is None:
            info = ''
        if type(info) is not type(' '):
            raise error('wrong info argument for MessageRedirector constructor')
        elif info in {'stdout', 'stderr', 'cout', 'cerr'}:
            self.info = pystir.newTextPrinter(info)
            self.info_case = 0
        else:
            self.info = pystir.newTextWriter(info)
            self.info_case = 1
        pystir.openChannel(0, self.info)

        if warn is None:
            warn = ''
        if type(warn) is not type(' '):
            raise error('wrong warn argument for MessageRedirector constructor')
        elif warn in {'stdout', 'stderr', 'cout', 'cerr'}:
            self.warn = pystir.newTextPrinter(warn)
            self.warn_case = 0
        else:
            self.warn = pystir.newTextWriter(warn)
            self.warn_case = 1
        pystir.openChannel(1, self.warn)

        if errr is None:
            errr = ''
        if type(errr) is not type(' '):
            raise error('wrong errr argument for MessageRedirector constructor')
        elif errr in {'stdout', 'stderr', 'cout', 'cerr'}:
            self.errr = pystir.newTextPrinter(errr)
            self.errr_case = 0
        else:
            self.errr = pystir.newTextWriter(errr)
            self.errr_case = 1
        pystir.openChannel(2, self.errr)

    def __del__(self):
        if self.info_case == 0:
            try_calling(pystir.deleteTextPrinter(self.info))
        else:
            try_calling(pystir.deleteTextWriter(self.info))
        pystir.closeChannel(0, self.info)
        if self.warn_case == 0:
            try_calling(pystir.deleteTextPrinter(self.warn))
        else:
            try_calling(pystir.deleteTextWriter(self.warn))
        pystir.closeChannel(1, self.warn)
        if self.errr_case == 0:
            try_calling(pystir.deleteTextPrinter(self.errr))
        else:
            try_calling(pystir.deleteTextWriter(self.errr))
        pystir.closeChannel(2, self.errr)

class Shape:
    '''
    Class for an abstract geometric shape used as a building block for
    creating phantom images.
    '''
    def __init__(self):
        self.handle = None
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def set_origin(self, origin):
        '''
        Sets the (discrete) coordinates of the shape centre on a voxel grid.
        '''
        _set_float_par(self.handle, 'Shape', 'x', origin[0])
        _set_float_par(self.handle, 'Shape', 'y', origin[1])
        _set_float_par(self.handle, 'Shape', 'z', origin[2])
    def get_origin(self):
        '''
        Returns the coordinates of the shape centre on a voxel grid.
        '''
        x = _float_par(self.handle, 'Shape', 'x')
        y = _float_par(self.handle, 'Shape', 'y')
        z = _float_par(self.handle, 'Shape', 'z')
        return (x, y, z)

class EllipticCylinder(Shape):
    '''
    Class for elliptic cylinder shape.
    '''
    def __init__(self):
        self.handle = None
        self.name = 'EllipsoidalCylinder'
        self.handle = pystir.cSTIR_newObject(self.name)
        check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def set_length(self, value):
        _set_float_par(self.handle, self.name, 'length', value)
    def get_length(self):
        return _float_par(self.handle, self.name, 'length')
    def set_radius_x(self, value):
        _set_float_par(self.handle, self.name, 'radius_x', value)
    def get_radius_x(self):
        return _float_par(self.handle, self.name, 'radius_x')
    def set_radius_y(self, value):
        _set_float_par(self.handle, self.name, 'radius_y', value)
    def get_radius_y(self):
        return _float_par(self.handle, self.name, 'radius_y')
    def set_radii(self, radii):
        _set_float_par(self.handle, self.name, 'radius_x', radii[0])
        _set_float_par(self.handle, self.name, 'radius_y', radii[1])
    def get_radii(self):
        rx = _float_par(self.handle, self.name, 'radius_x')
        ry = _float_par(self.handle, self.name, 'radius_y')
        return (rx, ry)

class DataContainer(ABC):
    '''
    Abstract base class for an abstract data container.
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

        Since this class is abstract, its methods cannot itself create a new
        object when e.g. adding two objects of this class, so new object is
        created by the first object using its same_object() method - see
        __add__ below.
        '''
        pass
    def norm(self):
        '''
        Returns the 2-norm of the container data viewed as a vector.
        '''
        assert self.handle is not None
        handle = pystir.cSTIR_norm(self.handle)
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
        assert_validities(self, other)
        handle = pystir.cSTIR_dot(self.handle, other.handle)
        check_status(handle)
        r = pyiutil.floatDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return r
    def __add__(self, other):
        '''
        Overloads + for data containers.

        Returns the sum of the container data with another container 
        data viewed as vectors.
        other: DataContainer
        '''
        assert_validities(self, other)
        z = self.same_object()
        z.handle = pystir.cSTIR_axpby(1.0, self.handle, 1.0, other.handle)
        check_status(z.handle)
        return z;
    def __sub__(self, other):
        '''
        Overloads - for data containers.

        Returns the difference of the container data with another container 
        data viewed as vectors.
        other: DataContainer
        '''
        assert_validities(self, other)
        z = self.same_object()
        z.handle = pystir.cSTIR_axpby(1.0, self.handle, -1.0, other.handle)
        check_status(z.handle)
        return z;
    def __mul__(self, other):
        '''
        Overloads * for data containers multiplication by a scalar or another
        data container.

        Returns the product self*other if other is a scalar
        or the dot product if it is DataContainer.
        other: DataContainer or a (real or complex) scalar
        '''
        assert self.handle is not None
        if type(self) == type(other):
            return self.dot(other)
        z = self.same_object()
        if type(other) == type(0):
            other = float(other)
        if type(other) == type(0.0):
            z.handle = pystir.cSTIR_mult(other, self.handle)
            z.src = 'mult'
            check_status(z.handle)
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
        if type(other) == type(0.0):
            z.handle = pystir.cSTIR_mult(other, self.handle)
            check_status(z.handle)
            return z;
        else:
            raise error('wrong multiplier')

class ImageData(DataContainer):
    '''Class for PET image data objects.

    ImageData objects contains both geometric data and the actual voxel
    values. You have to use the `as_array` method to get an array with
    the voxel values, and use the `fill` function to change the voxel values.
    '''
    def __init__(self, arg = None):
        '''Create an ImageData object

        Arguments:
        str            : read the object from a file specified by <arg>
                         (the file format has to be support by STIR).
        AcquisitionData: create an object compatible with the scanner data
                         recorded in an AcquisitionData object <arg>.
                         This sets default voxel sizes.
        None           : create an empty ImageData object. Call initialise()
                         method before using it.
        '''
        self.handle = None
        if isinstance(arg, str):
            self.handle = pystir.cSTIR_objectFromFile('Image', arg)
            check_status(self.handle)
        elif isinstance(arg, AcquisitionData):
            assert arg.handle is not None
            self.handle = pystir.cSTIR_imageFromAcquisitionData(arg.handle)
            check_status(self.handle)
        elif arg is not None:
            raise error\
                  ('wrong argument ' + repr(arg) + ' for ImageData constructor')
        self.name = 'ImageData'
        self.rimsize = -1
    def __del__(self):
        '''Deallocates this ImageData object.'''
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def same_object(self):
        '''See DataContainer.same_object().
        '''
        return ImageData()
    def initialise\
        (self, arg1, arg2 = 0, arg3 = 0, arg4 = 1, arg5 = 1, arg6 = 1, \
         arg7 = 0, arg8 = 0, arg9 = 0):
        '''Change image size and geometric information

        Sets this image size in voxels, voxel sizes in mm and the origin.
        All arguments except the first one are optional.
        Present arguments are either all scalars or all tuples.
        The first tuple argument or three scalar arguments set the image
        sizes in voxels.
        The second tuple argument or three scalar arguments set the voxel
        sizes in mm (if absent, sizes default to (1,1,1)).
        The third tuple argument or three scalar arguments set the origin
        (if absent, defaults to (0,0,0)).
        '''
        if type(arg1) == type((0,0,0)):
            dim = arg1
            if arg2 == 0:
                vsize = (1, 1, 1)
                origin = (0, 0, 0)
            else:
                vsize = arg2
                if arg3 == 0:
                    origin = (0.0, 0.0, 0.0)
                else:
                    origin = arg3
        else:
            dim = (arg1, arg2, arg3)
            vsize = (arg4, arg5, arg6)
            origin = (arg7, arg8, arg9)
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
        self.handle = None
        voxels = pystir.cSTIR_voxels3DF\
                      (dim[0], dim[1], dim[2], \
                       vsize[0], vsize[1], vsize[2], \
                       origin[0], origin[1], origin[2])
        check_status(voxels)
        self.handle = pystir.cSTIR_imageFromVoxels(voxels)
        check_status(self.handle)
        pyiutil.deleteDataHandle(voxels)
    def fill(self, value):
        '''Sets the voxel-values.

        The argument is either 3D Numpy ndarray of values or a scalar to be
        assigned at each voxel. When using an ndarray, the array size has to
        have the same size as an array returned by `as_array`.
        '''
        assert self.handle is not None
        if isinstance(value, numpy.ndarray):
            if value.dtype is numpy.dtype('float32'):
                #print('keeping dtype float32')
                v = value
            else:
                #print('changing dtype to float32')
                v = value.astype(numpy.float32)
            try_calling(pystir.cSTIR_setImageData(self.handle, v.ctypes.data))
        elif isinstance(value, float):
            try_calling(pystir.cSTIR_fillImage(self.handle, value))
        elif isinstance(value, int):
            try_calling(pystir.cSTIR_fillImage(self.handle, float(value)))
        else:
            raise error('wrong fill value.' + \
                        ' Should be numpy.ndarray, float or int')
        return self
    def clone(self):
        '''Creates a copy of this image.'''
        assert self.handle is not None
        image = ImageData()
        image.handle = pystir.cSTIR_imageFromImage(self.handle)
        check_status(image.handle)
        return image
    def get_uniform_copy(self, value = 1.0):
        '''Creates a copy of this image filled with <value>.'''
        assert self.handle is not None
        image = ImageData()
        image.handle = pystir.cSTIR_imageFromImage(self.handle)
        check_status(image.handle)
        image.fill(value)
        return image
    def add_shape(self, shape, scale):
        '''Adds a shape to self - see Shape above.
        '''
        assert self.handle is not None
        assert_validity(shape, Shape)
        try_calling(pystir.cSTIR_addShape(self.handle, shape.handle, scale))
    def read_from_file(self, filename):
        '''Read data from file.

        Replaces the current content of the object.
        '''
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
        self.handle = pystir.cSTIR_objectFromFile('Image', filename)
        check_status(self.handle)
    def write(self, filename):
        '''Writes self to an Interfile - see STIR documentation for details.
        '''
        assert self.handle is not None
        try_calling(pystir.cSTIR_writeImage(self.handle, filename))
    def dimensions(self):
        '''Returns image dimensions as a tuple (nx, ny, nz).'''
        assert self.handle is not None
        dim = numpy.ndarray((3,), dtype = numpy.int32)
        try_calling \
            (pystir.cSTIR_getImageDimensions(self.handle, dim.ctypes.data))
        return tuple(dim[::-1])
    def voxel_sizes(self):
        '''Returns image voxel sizes as a tuple (vx, vy, vz).'''
        assert self.handle is not None
        vs = numpy.ndarray((3,), dtype = numpy.float32)
        try_calling \
            (pystir.cSTIR_getImageVoxelSizes(self.handle, vs.ctypes.data))
        return tuple(vs[::-1])
    def as_array(self):
        '''Returns 3D Numpy ndarray with values at the voxels.'''
        assert self.handle is not None
        dim = numpy.ndarray((9,), dtype = numpy.int32)
        try_calling \
            (pystir.cSTIR_getImageDimensions(self.handle, dim.ctypes.data))
        nz = dim[0]
        ny = dim[1]
        nx = dim[2]
        if nx == 0 or ny == 0 or nz == 0:
            raise error('image data not available')
        array = numpy.ndarray((nz, ny, nx), dtype = numpy.float32)
        try_calling(pystir.cSTIR_getImageData(self.handle, array.ctypes.data))
        return array
    def show(self):
        '''Displays xy-cross-sections of this image at z selected interactively.'''
        assert self.handle is not None
        if not HAVE_PYLAB:
            print('pylab not found')
            return
        data = self.as_array()
        nz = data.shape[0]
        print('Please enter slice numbers (e.g.: 1, 3-5)')
        print('(a value outside the range [1 : %d] will stop this loop)' % nz)
        while True:
            s = str(input('slices to display: '))
            if len(s) < 1:
                break
            err = show_3D_array(data, index = s, label = 'slice')
            if err != 0:
                print('out-of-range slice numbers selected, quitting the loop')
                break

DataContainer.register(ImageData)

class ImageDataProcessor:
    '''Class for image processors.

    An ImageDataProcessor changes an image in some way, e.g. by filtering.'''
    def __init__(self):
        self.handle = None
        # TODO: handle input and output in cSTIR
        self.input = None
        self.output = None
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def apply(self, image):
        '''Applies this filter to the specified image.'''
        assert_validity(image, ImageData)
        try_calling(pystir.cSTIR_applyImageDataProcessor\
                 (self.handle, image.handle))
    def set_input(self, input):
        '''
        Sets the input data.
        '''
        assert isinstance(input, ImageData)
        self.input = input
    def process(self, input = None):
        '''
        Processes data.
        '''
        if input is not None:
            self.input = input
        if self.input is None:
            raise error('input image not set')
        assert isinstance(self.input, ImageData)
        self.output = self.input.clone()
        self.apply(self.output)
        return self.output
    def get_output(self):
        '''
        Returns the output data.
        '''
        return self.output

class TruncateToCylinderProcessor(ImageDataProcessor):
    '''
    Class for the image filter that zeroes the image outside the cylinder
    of the same xy-diameter and z-size as those of the image.
    '''
    def __init__(self):
        self.handle = None
        self.name = 'TruncateToCylindricalFOVImageProcessor'
        self.handle = pystir.cSTIR_newObject(self.name)
        check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def set_strictly_less_than_radius(self, flag):
        '''
        Specifies whether the area not affected by filtering is strictly
        inside the cylinder (flag = True) or not (flag = False).
        '''
        _set_char_par\
            (self.handle, 'TruncateToCylindricalFOVImageProcessor',\
             'strictly_less_than_radius', repr(flag))
    def get_strictly_less_than_radius(self):
        '''
        Returns the answer to the question: Is the area not affected by 
        filtering strictly inside the cylinder?
        '''
        return _int_par\
               (self.handle, 'TruncateToCylindricalFOVImageProcessor',\
                'strictly_less_than_radius') != 0

class RayTracingMatrix:
    '''
    Class for objects holding sparse matrix representation of the ray
    tracing projector G (see AcquisitionModel class).
    '''
    def __init__(self):
        self.handle = None
        self.name = 'RayTracingMatrix'
        self.handle = pystir.cSTIR_newObject(self.name)
        check_status(self.handle)
        _set_int_par(self.handle, self.name, 'num_tangential_LORs', 2)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def set_num_tangential_LORs(self, value):
        '''
        Set the number of LORs (or rays) for each bin in the sinogram.
        They are currently (approximately) parallel and spaced in the
        tangential direction (i.e. orthogonal to the axial direction).
        '''
        _set_int_par(self.handle, self.name, 'num_tangential_LORs', value)
        return self
    def get_num_tangential_LORs(self):
        '''
        Returns the number of LORs for each bin in the sinogram.
        '''
        return _int_par(self.handle, self.name, 'num_tangential_LORs')

class AcquisitionData(DataContainer):
    '''Class for PET acquisition data.'''
    def __init__\
        (self, src = None, span = 1, max_ring_diff = -1, view_mash_factor = 1):
        ''' 
        Creates new AcquisitionData object from a file or another
        AcquisitionData object;
        src:  file name (Python str) or AcquisitionData object or scanner name
        '''
        self.handle = None
        self.name = 'AcquisitionData'
        self.read_only = False
        self.src = None
        if src is None:
            return
        if isinstance(src, str):
            i = src.find('.')
            if i > -1:
                # src is a file name
                self.handle = pystir.cSTIR_objectFromFile\
                            ('AcquisitionData', src)
                self.read_only = True
                self.src = 'file'
            else:
                # src is a scanner name
                self.handle = pystir.cSTIR_acquisitionsDataFromScannerInfo\
                    (src, span, max_ring_diff, view_mash_factor)
                if pyiutil.executionStatus(self.handle) != 0:
                    msg = pyiutil.executionError(self.handle)
                    if msg == 'Unknown scanner':
                        raise error\
                            ('Unknown scanner ' + src + \
                             ' or missing raw data file extension')
                self.src = 'scanner'
        elif isinstance(src, AcquisitionData):
            # src is AcquisitionData
            assert src.handle is not None
            self.handle = pystir.cSTIR_acquisitionsDataFromTemplate\
                (src.handle)
            self.src = 'template'
        else:
            raise error('Wrong source in AcquisitionData constructor')
        check_status(self.handle)
    def __del__(self):
        #print('deleting AcquisitionData object originated from ', self.src)
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
        try_calling(pystir.cSTIR_setAcquisitionsStorageScheme(scheme))
    @staticmethod
    def get_storage_scheme():
        '''Returns acquisition data storage scheme.
        '''
        handle = pystir.cSTIR_getAcquisitionsStorageScheme()
        check_status(handle)
        scheme = pyiutil.charDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return scheme
    def same_object(self):
        '''See DataContainer.same_object().
        '''
        return AcquisitionData()
    def read_from_file(self, filename): # 'read_from_file' is misleading
        '''Read data from file.

        Replaces the current content of the object.
        '''
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
        self.handle = pystir.cSTIR_objectFromFile('AcquisitionData', filename)
        check_status(self.handle)
        self.read_only = True
    def create_uniform_image(self, value = 0, xy = None):
        ''' 
        Creates ImageData object containing PET image of dimensions
        and voxel sizes compatible with the scanner geometry stored
        in this AcquisitionData object and assigns a given value
        to all voxels;
        value:  a Python float.
        xy: y and x dimensions tuple
        '''
        assert self.handle is not None
        image = ImageData()
        if xy is None:
            image.handle = pystir.cSTIR_imageFromAcquisitionData(self.handle)
        else:
            image.handle = pystir.cSTIR_imageFromAcquisitionDataAndNxNy\
                           (self.handle, xy[1], xy[0])
        check_status(image.handle)
        image.fill(value)
        return image
    def dimensions(self):
        ''' Returns a tuple of the data dimensions:
        - number of sinograms
        - number of views
        - number of tangential positions.
        '''
        assert self.handle is not None
        dim = numpy.ndarray((3,), dtype = numpy.int32)
        try_calling(pystir.cSTIR_getAcquisitionsDimensions\
            (self.handle, dim.ctypes.data))
        nt = dim[0]
        nv = dim[1]
        ns = dim[2]
        return ns, nv, nt
    def as_array(self):
        ''' 
        Returns a copy of acquisition data stored in this object as a
        NumPy ndarray of 3 dimensions (in default C ordering of data):
        - number of sinograms
        - number of views
        - number of tangential positions.
        '''
        assert self.handle is not None
        dim = numpy.ndarray((3,), dtype = numpy.int32)
        try_calling(pystir.cSTIR_getAcquisitionsDimensions\
            (self.handle, dim.ctypes.data))
        nx = dim[0]
        ny = dim[1]
        nz = dim[2]
        if nx == 0 or ny == 0 or nz == 0:
            raise error('Acquisition data not present')
        array = numpy.ndarray((nz, ny, nx), dtype = numpy.float32)
        try_calling(pystir.cSTIR_getAcquisitionsData\
            (self.handle, array.ctypes.data))
        return array
    def fill(self, value):
        ''' 
        Fills the object with values;
        value:  either NumPy ndarray or another AcquisitionData object
                or Python float.
        '''
        assert self.handle is not None
        if self.read_only:
            raise error('Cannot fill read-only object, consider filling a clone')
        if isinstance(value, numpy.ndarray):
            if value.dtype is numpy.dtype('float32'):
                #print('keeping dtype float32')
                v = value
            else:
                #print('changing dtype to float32')
                v = value.astype(numpy.float32)
            try_calling(pystir.cSTIR_setAcquisitionsData\
                        (self.handle, v.ctypes.data))
        elif isinstance(value, AcquisitionData):
            assert value.handle is not None
            try_calling(pystir.cSTIR_fillAcquisitionsDataFromAcquisitionsData\
                (self.handle, value.handle))
        elif isinstance(value, float):
            try_calling(pystir.cSTIR_fillAcquisitionsData(self.handle, value))
        elif isinstance(value, int):
            try_calling(pystir.cSTIR_fillAcquisitionsData\
                        (self.handle, float(value)))
        else:
            raise error('Wrong fill value.' + \
                ' Should be numpy.ndarray, AcquisitionData, float or int')
        return self
    def write(self, filename):
        '''Writes self to an Interfile - see STIR documentation for details.
        '''
        assert self.handle is not None
        try_calling(pystir.cSTIR_writeAcquisitionData(self.handle, filename))
    def clone(self):
        ''' 
        Returns a true copy of this object (not Python handle).
        '''
        ad = AcquisitionData(self)
        ad.fill(self)
        ad.src = 'clone'
        return ad
    def get_uniform_copy(self, value = 0):
        ''' 
        Returns a true copy of this object filled with a given value;
        value:  a Python float.
        '''
        ad = AcquisitionData(self)
        ad.fill(value)
        ad.src = 'copy'
        return ad
    def rebin(self, num_segments_to_combine, \
        num_views_to_combine = 1, num_tang_poss_to_trim = 0, \
        do_normalisation = True, max_in_segment_num_to_process = -1):
        ad = AcquisitionData()
        ad.handle = pystir.cSTIR_rebinnedAcquisitionData(self.handle, \
            num_segments_to_combine, num_views_to_combine, \
            num_tang_poss_to_trim, do_normalisation, \
            max_in_segment_num_to_process)
        check_status(ad.handle)
        return ad

DataContainer.register(AcquisitionData)

class ListmodeToSinograms:
    '''
    Class for listmode-to-sinogram converter.

    This class reads list mode data and produces corresponding *sinograms*,
    i.e. histogrammed data in the format of PETAcquisitionData.

    It has two main functions:
      - process() can be used to read prompts and/or delayed coincidences to
        produce a single PETAcquisitionData.
        Two conversion flags decide what is to be done with 3 possible cases:
        - `store_prompts`=`true`, `store_delayeds`=`false`: only prompts stored
        - `store_prompts`=`false`, `store_delayeds`=`true`: only delayeds stored
        - `store_prompts`=`true`, `store_delayeds`=`true`: prompts-delayeds stored
        Clearly, enabling the `store_delayeds` option only makes sense if the
        data was acquired accordingly.
      - estimate_randoms() can be used to get a relatively noiseless estimate of the 
        random coincidences.

    Currently, the randoms are estimated from the delayed coincidences using the
    following strategy:
       1. singles (one per detector) are estimated using a Maximum Likelihood
          estimator
       2. randoms-from-singles are computed per detector-pair via the usual
          product formula. These are then added together for all detector pairs
          in a certain histogram-bin in the data (accommodating for view mashing
          and axial compression).

    The actual algorithm is described in

    D. Hogg, K. Thielemans, S. Mustafovic, and T. J. Spinks,
    "A study of bias for various iterative reconstruction methods in PET,"
    in 2002 IEEE Nuclear Science Symposium Conference Record, vol. 3. IEEE,
    Nov. 2002, pp. 1519-1523 (http://dx.doi.org/10.1109/nssmic.2002.1239610).
    '''
    def __init__(self, file = None):
        self.handle = None
        self.name = 'ListmodeToSinograms'
        if file is None:
            self.handle = pystir.cSTIR_newObject(self.name)
        else:
            self.handle = pystir.cSTIR_objectFromFile(self.name, file)
        self.output = None
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def set_input(self, lm_file):
        '''Sets the listmode file name.
        '''
        _set_char_par(self.handle, self.name, 'input', lm_file)
    def set_output_prefix(self, sino_file):
        '''Sets the sinograms file names prefix.
        '''
        _set_char_par(self.handle, self.name, 'output', sino_file)
    def set_template(self, templ):
        '''Sets the sinograms template.
        '''
        _set_char_par(self.handle, self.name, 'template', templ)
    def set_time_interval(self, start, stop):
        '''Sets time interval.

        Only data scanned during this time interval will be converted.
        '''
        interval = numpy.ndarray((2,), dtype = numpy.float32)
        interval[0] = start
        interval[1] = stop
        try_calling(pystir.cSTIR_setListmodeToSinogramsInterval\
            (self.handle, interval.ctypes.data))
    def flag_on(self, flag):
        '''Switches on (sets to 'true') a conversion flag (see conversion flags
           description above).
        '''
        try_calling(pystir.cSTIR_setListmodeToSinogramsFlag\
            (self.handle, flag, 1))
    def flag_off(self, flag):
        '''Switches off (sets to 'false') a conversion flag (see conversion flags
           description above).
        '''
        try_calling(pystir.cSTIR_setListmodeToSinogramsFlag\
            (self.handle, flag, 0))
    def set_up(self):
        '''Sets up the conversion.
        '''
        try_calling(pystir.cSTIR_setupListmodeToSinogramsConverter(self.handle))
    def process(self):
        '''Performs the conversion.
        '''
        self.output = AcquisitionData()
        self.output.handle = \
                           pystir.cSTIR_convertListmodeToSinograms(self.handle)
        check_status(self.output.handle)
    def get_output(self):
        '''Returns the sinograms as an AcquisitionData object.
        '''
        if self.output is None:
            raise error('Conversion to sinograms not done')
        return self.output
    def estimate_randoms(self):
        '''Returns an estimate of the randoms as an AcquisitionData object.
        '''
        randoms = AcquisitionData()
        randoms.handle = pystir.cSTIR_computeRandoms(self.handle)
        check_status(randoms.handle)
        return randoms

class AcquisitionSensitivityModel:
    '''
    Class that handles PET scanner detector efficiencies and attenuation.

    Is used by AcquisitionModel (see below) for multiplication by 1/n.
    '''
    def __init__(self, src, other_src = None):
        ''' 
        Creates new AcquisitionSensitivityModel object
        - from an ECAT8 file or
        - from ImageData object containing attenuation image (units: 1/cm) or
        - from AcquisitionData object containing bin efficiencies or
        - by chaining two existing AcquisitionSensitivityModel objects
        src: file name or ImageData object or AcquisitionData object
        pther_src: AcquisitionSensitivityModel object (optional)
        '''
        self.handle = None
        self.name = 'AcquisitionSensitivityModel'
        if src is None:
            return
        if isinstance(src, str):
            # create from ECAT8 file
            handle = pyiutil.charDataHandle(src)
            self.handle = pystir.cSTIR_createPETAcquisitionSensitivityModel\
                (handle, 'n')
        elif isinstance(src, ImageData):
            # create from attenuation image
            assert src.handle is not None
            assert isinstance(other_src, AcquisitionModel)
            self.handle = pystir.cSTIR_createPETAttenuationModel\
                          (src.handle, other_src.handle)
        elif isinstance(src, AcquisitionData):
            # create from bin efficiencies (cf. AcquisitionModel)
            assert src.handle is not None
            self.handle = pystir.cSTIR_createPETAcquisitionSensitivityModel\
                (src.handle, 's')
        elif isinstance(src, AcquisitionSensitivityModel) and \
             isinstance(other_src, AcquisitionSensitivityModel):
            # chain two acquisition sensitivity models 
            assert src.handle is not None
            assert other_src.handle is not None
            self.handle = pystir.cSTIR_chainPETAcquisitionSensitivityModels\
                (src.handle, other_src.handle)
        else:
            raise error\
                  ('Wrong source in AcquisitionSensitivityModel constructor')
        check_status(self.handle)
    def set_up(self, ad):
        '''Sets up the object.
        '''
        assert self.handle is not None
        assert_validity(ad, AcquisitionData)
        try_calling(pystir.cSTIR_setupAcquisitionSensitivityModel\
            (self.handle, ad.handle))
    def normalise(self, ad):
        '''Multiplies the argument by n (cf. AcquisitionModel).
           If self is a chain of two AcquisitionSensitivityModels, then n is
           a product of two normalisations.
        '''
        assert self.handle is not None
        assert_validity(ad, AcquisitionData)
        try_calling(pystir.cSTIR_applyAcquisitionSensitivityModel\
            (self.handle, ad.handle, 'normalise'))
    def unnormalise(self, ad):
        '''Multiplies the argument by 1/n (cf. AcquisitionModel).
           If self is a chain of two AcquisitionSensitivityModels, then n is
           a product of two normalisations.
        '''
        assert self.handle is not None
        assert_validity(ad, AcquisitionData)
        try_calling(pystir.cSTIR_applyAcquisitionSensitivityModel\
            (self.handle, ad.handle, 'unnormalise'))
    def forward(self, ad):
        '''Same as unnormalise except that the argument remains unchanged
           and  a new AcquisitionData equal to the argument multiplied
           by 1/n is returned.
        '''
        assert self.handle is not None
        assert_validity(ad, AcquisitionData)
        fd = AcquisitionData()
        fd.handle = pystir.cSTIR_applyAcquisitionSensitivityModel\
            (self.handle, ad.handle, 'fwd')
        check_status(fd.handle)
        return fd
    def invert(self, ad):
        '''Same as normalise except that the argument remains unchanged
           and  a new AcquisitionData equal to the argument multiplied
           by n is returned.
        '''
        assert self.handle is not None
        assert_validity(ad, AcquisitionData)
        fd = AcquisitionData()
        fd.handle = pystir.cSTIR_applyAcquisitionSensitivityModel\
            (self.handle, ad.handle, 'inv')
        check_status(fd.handle)
        return fd
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

class AcquisitionModel:
    ''' 
    Class for a PET acquisition model that relates an image x to the
    acquisition data y as
    (F)    y = [1/n] (G x + [a]) + [b]
    where:
    G is the geometric (ray tracing) projector from the image voxels
    to the scanner's pairs of detectors (bins);
    a and b are otional additive and background terms representing
    the effects of noise and scattering; assumed to be 0 if not present;
    n is an optional bin normalization term representing the inverse of
    detector (bin) efficiencies; assumed to be 1 if not present.
    The computation of y for a given x by the above formula (F) is
    referred to as forward projection, and the computation of
    (B)    z = G' m y
    where G' is the transpose of G and m = 1/n, is referred to as 
    backward projection.
    '''
    def __init__(self):
        self.handle = None
        self.name = 'AcquisitionModel'
    def set_up(self, acq_templ, img_templ):
        ''' 
        Prepares this object for performing forward and backward
        projections;
        acq_templ:  an AcquisitionData object used as a template for
                    creating an AcquisitionData object to store forward
                    projection;
        img_templ:  an ImageData object used as a template for creating an
                    ImageData object to store backward projection.
        '''
        assert_validity(acq_templ, AcquisitionData)
        assert_validity(img_templ, ImageData)
        try_calling(pystir.cSTIR_setupAcquisitionModel\
            (self.handle, acq_templ.handle, img_templ.handle))
    def set_additive_term(self, at):
        ''' 
        Sets the additive term a in (F);
        at:  an AcquisitionData object containing a.
        '''
        assert_validity(at, AcquisitionData)
        _setParameter\
            (self.handle, 'AcquisitionModel', 'additive_term', at.handle)
    def set_background_term(self, bt):
        ''' 
        Sets the background term b in (F);
        bt:  an AcquisitionData object containing b.
        '''
        assert_validity(bt, AcquisitionData)
        _setParameter\
            (self.handle, 'AcquisitionModel', 'background_term', bt.handle)
    def set_acquisition_sensitivity(self, asm):
        ''' 
        Sets the normalization n in (F);
        norm:  an AcquisitionSensitivityModel object containing normalisation n
        '''
        assert_validity(asm, AcquisitionSensitivityModel)
        _setParameter\
            (self.handle, 'AcquisitionModel', 'asm', asm.handle)
    def forward(self, image):
        ''' 
        Returns the forward projection of x given by (F);
        image   :  an ImageData object containing x;
        '''
        assert_validity(image, ImageData)
        ad = AcquisitionData()
        ad.handle = pystir.cSTIR_acquisitionModelFwd(self.handle, image.handle)
        check_status(ad.handle)
        return ad;
    def backward(self, ad):
        ''' 
        Returns the backward projection of y giben by (B);
        ad:  an AcquisitionData object containing y.
        '''
        assert_validity(ad, AcquisitionData)
        image = ImageData()
        image.handle = pystir.cSTIR_acquisitionModelBwd\
            (self.handle, ad.handle)
        check_status(image.handle)
        return image

class AcquisitionModelUsingMatrix(AcquisitionModel):
    ''' 
    Class for a PET acquisition model that uses (implicitly) a sparse
    matrix for G in (F).
    '''
    def __init__(self, matrix = None):
        ''' 
        Creates an AcquisitionModelUsingMatrix object, optionally setting
        the ray tracing matrix to be used for projecting;
        matrix:  a RayTracingMatrix object to represent G in (F).
        '''
        self.handle = None
        self.name = 'AcqModUsingMatrix'
        self.handle = pystir.cSTIR_newObject(self.name)
        check_status(self.handle)
        if matrix is None:
            matrix = RayTracingMatrix()
        assert_validity(matrix, RayTracingMatrix)
        _setParameter(self.handle, self.name, 'matrix', matrix.handle)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def set_matrix(self, matrix):
        ''' 
        Sets the ray tracing matrix to be used for projecting;
        matrix:  a RayTracingMatrix object to represent G in (F).
        '''
        assert_validity(matrix, RayTracingMatrix)
        _setParameter(self.handle, self.name, 'matrix', matrix.handle)
##    def get_matrix(self):
##        ''' 
##        Returns the ray tracing matrix used for projecting;
##        matrix:  a RayTracingMatrix object representing G in (F).
##        '''
##        matrix = RayTracingMatrix()
##        matrix.handle = pystir.cSTIR_parameter(self.handle, self.name, 'matrix')
##        check_status(matrix.handle)
##        return matrix

class AcquisitionModelUsingRayTracingMatrix(AcquisitionModelUsingMatrix):
    ''' 
    Class for a PET acquisition model that uses (implicitly) a ray tracing
    matrix for G in (F).
    '''
    def __init__(self, matrix = None):
        ''' 
        Creates an AcquisitionModelUsingMatrix object, optionally setting
        the ray tracing matrix to be used for projecting;
        matrix:  a RayTracingMatrix object to represent G in (F).
        '''
        self.handle = None
        self.name = 'AcqModUsingMatrix'
        self.handle = pystir.cSTIR_newObject(self.name)
        check_status(self.handle)
        if matrix is None:
            matrix = RayTracingMatrix()
        assert_validity(matrix, RayTracingMatrix)
        self.matrix = matrix
        _setParameter(self.handle, self.name, 'matrix', matrix.handle)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
##    def set_matrix(self, matrix):
##        ''' 
##        Sets the ray tracing matrix to be used for projecting;
##        matrix:  a RayTracingMatrix object to represent G in (F).
##        '''
##        _setParameter(self.handle, self.name, 'matrix', matrix.handle)
    def get_matrix(self):
        ''' 
        Returns the ray tracing matrix used for projecting;
        matrix:  a RayTracingMatrix object representing G in (F).
        '''
        matrix = RayTracingMatrix()
        matrix.handle = pystir.cSTIR_parameter(self.handle, self.name, 'matrix')
        check_status(matrix.handle)
        return matrix
    def set_num_tangential_LORs(self, value):
        '''
        Set the number of LORs (or rays) for each bin in the sinogram.
        They are currently (approximately) parallel and spaced in the
        tangential direction (i.e. orthogonal to the axial direction).
        '''
##        return self.matrix.set_num_tangential_LORs(value)
        return self.get_matrix().set_num_tangential_LORs(value)

class Prior:
    '''
    Class for objects handling the prior: a penalty term to be added to the
    objective function maximized by iterative reconstruction algorithms.
    '''
    def __init__(self):
        self.handle = None
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def set_penalisation_factor(self, value):
        '''
        Sets the factor by which the penalty term (prior) is to be multiplied
        before adding to the objective function.
        '''
        _set_float_par\
            (self.handle, 'GeneralisedPrior', 'penalisation_factor', value)
        return self
    def get_penalisation_factor(self):
        '''
        Returns the penalty factor in front of the prior.
        '''
        return _float_par\
            (self.handle, 'GeneralisedPrior', 'penalisation_factor')
    def get_gradient(self, image):
        '''
        Returns the value of the gradient of the prior for the specified image.
        image: ImageData object
        '''
        assert isinstance(image, ImageData)
        grad = ImageData()
        grad.handle = pystir.cSTIR_priorGradient(self.handle, image.handle)
        check_status(grad.handle)
        return grad
    def set_up(self):
        try_calling(pystir.cSTIR_setupPrior(self.handle))
##        handle = pystir.cSTIR_setupObject('GeneralisedPrior', self.handle)
##        check_status(handle)
##        pyiutil.deleteDataHandle(handle)

class QuadraticPrior(Prior):
    '''
    Class for the prior that is a quadratic functions of the image values.
    '''
    def __init__(self):
        self.handle = None
        self.name = 'QuadraticPrior'
        self.handle = pystir.cSTIR_newObject(self.name)
        check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

class PLSPrior(Prior):
    '''
    Class for PLS prior.
    '''
    def __init__(self):
        self.handle = None
        self.name = 'PLSPrior'
        self.handle = pystir.cSTIR_newObject(self.name)
        check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def set_anatomical_image(self, image):
        assert isinstance(image, ImageData)
        _setParameter(self.handle, 'PLSPrior',\
            'anatomical_image', image.handle)

class ObjectiveFunction:
    '''
    Class for the objective function maximized by the iterative reconstruction
    algorithms.
    '''
    def __init__(self):
        self.handle = None
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def set_prior(self, prior):
        '''
        Sets the prior (penalty term to be added to the objective function).
        '''
        assert_validity(prior, Prior)
        _setParameter(self.handle, 'GeneralisedObjectiveFunction',\
            'prior', prior.handle)
        self.prior = prior
    def get_prior(self):
        '''
        Returns the prior currently used by this objective function.
        '''
        prior = Prior()
        prior.handle = pystir.cSTIR_parameter\
            (self.handle, 'GeneralisedObjectiveFunction', 'prior')
        check_status(prior.handle)
        return prior
    def set_num_subsets(self, n):
        '''
        Sets the number of subsets of ray projections to be used for computing
        additive components of the gradient used by Ordered Subset algorithms
        for maximizing this objective function.
        If the ray tracing projector G is a matrix, the subsets in question are 
        subsets of its rows.
        n: number of subsets, Python integer scalar
        '''
        _set_int_par\
            (self.handle, 'GeneralisedObjectiveFunction', 'num_subsets', n)
    def set_up(self, image):
        '''
        Prepares this object for use.
        image: ImageData object
        '''
        assert_validity(image, ImageData)
        try_calling(pystir.cSTIR_setupObjectiveFunction\
                    (self.handle, image.handle))
    def value(self, image):
        '''
        Returns the value of this objective function on the specified image.
        image: ImageData object
        '''
        assert_validity(image, ImageData)
        handle = pystir.cSTIR_objectiveFunctionValue(self.handle, image.handle)
        check_status(handle)
        v = pyiutil.floatDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return v
    def get_value(self, image):
        '''
        Returns the value of this objective function on the specified image.
        image: ImageData object
        '''
        return self.value(image)
    def gradient(self, image, subset = -1):
        '''
        Returns the value of the additive component of the gradient of this 
        objective function on the specified image corresponding to the specified
        subset (see set_num_subsets() method).
        If no subset is specified, returns the full gradient, i.e. the sum of
        the subset components.
        image: ImageData object
        subset: Python integer scalar
        '''
        assert_validity(image, ImageData)
        grad = ImageData()
        grad.handle = pystir.cSTIR_objectiveFunctionGradient\
            (self.handle, image.handle, subset)
        check_status(grad.handle)
        return grad
    def get_gradient(self, image):
        '''
        Returns the gradient of the objective function on specified image.
        image: ImageData object
        '''
        return self.gradient(image)
    def get_subset_gradient(self, image, subset):
        '''
        Returns the value of the additive component of the gradient of this 
        objective function on <image> corresponding to the specified subset
        (see set_num_subsets() method).
        image: ImageData object
        subset: Python integer scalar
        '''
        return self.gradient(image, subset)

class PoissonLogLikelihoodWithLinearModelForMean(ObjectiveFunction):
    '''
    Class for STIR PoissonLogLikelihoodWithLinearModelForMean object, see
    http://stir.sourceforge.net/documentation/doxy/html/classstir_1_1PoissonLogLikelihoodWithLinearModelForMean.html
    '''
    def __init__(self):
        self.handle = None
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
##    def set_sensitivity_filename(self, name):
##        _set_char_par\
##            (self.handle, 'PoissonLogLikelihoodWithLinearModelForMean',\
##             'sensitivity_filename', name)
##    def set_use_subset_sensitivities(self, flag):
##        _set_char_par\
##            (self.handle, 'PoissonLogLikelihoodWithLinearModelForMean',\
##             'use_subset_sensitivities', repr(flag))
    def set_recompute_sensitivity(self, flag):
        _set_char_par\
            (self.handle, 'PoissonLogLikelihoodWithLinearModelForMean',\
             'recompute_sensitivity', repr(flag))
    def get_subset_sensitivity(self, subset):
        '''
        Returns an ImageData object containing sensitivity image for the
        specified subset.
        '''
        ss = ImageData()
        ss.handle = pystir.cSTIR_subsetSensitivity(self.handle, subset)
        check_status(ss.handle)
        return ss
    def get_backprojection_of_acquisition_ratio(self, image, subset):
        '''
        Computes back-projection of the ratio of measured to estimated 
        acquisition data.
        '''
        assert_validity(image, ImageData)
        grad = ImageData()
        grad.handle = pystir.cSTIR_objectiveFunctionGradientNotDivided\
            (self.handle, image.handle, subset)
        check_status(grad.handle)
        return grad

class PoissonLogLikelihoodWithLinearModelForMeanAndProjData\
(PoissonLogLikelihoodWithLinearModelForMean):
    '''
    Class for STIR PoissonLogLikelihoodWithLinearModelForMeanAndProjData object,
    see
    http://stir.sourceforge.net/documentation/doxy/html/classstir_1_1PoissonLogLikelihoodWithLinearModelForMeanAndProjData.html
    '''
    def __init__(self):
        self.handle = None
        self.name = 'PoissonLogLikelihoodWithLinearModelForMeanAndProjData'
        self.handle = pystir.cSTIR_newObject(self.name)
        check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def set_input_filename(self, name):
        '''Sets the name of the file containing acquisition data.
        '''
        _set_char_par\
            (self.handle, self.name, 'input_filename', name)
##    def set_zero_seg0_end_planes(self, flag):
##        _set_char_par\
##            (self.handle, self.name, 'zero_seg0_end_planes', repr(flag))
##    def set_max_segment_num_to_process(self, n):
##        _set_int_par(self.handle, self.name, 'max_segment_num_to_process', n)
    def set_acquisition_model(self, am):
        '''
        Sets the acquisition model to be used by this objective function.
        '''
        assert_validity(am, AcquisitionModel)
        _setParameter\
            (self.handle, self.name, 'acquisition_model', am.handle)
##    def get_acquisition_model(self):
##        '''
##        Returns the acquisition model used by this objective function.
##        '''
##        am = AcquisitionModelUsingMatrix()
##        if am.handle is not None:
##            pyiutil.deleteDataHandle(am.handle)
##        am.handle = pystir.cSTIR_parameter\
##            (self.handle, self.name, 'acquisition_model')
##        check_status(am.handle)
##        return am
    def set_acquisition_data(self, ad):
        '''
        Sets the acquisition data to be used by this objective function.
        '''
        assert_validity(ad, AcquisitionData)
        _setParameter\
            (self.handle, self.name, 'acquisition_data', ad.handle)

class Reconstructor:
    '''
    Class for a generic PET reconstructor.
    '''
    def __init__(self):
        self.handle = None
        self.input = None
        self.image = None
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
##    def set_input(self, input_data):
##        self.input = input_data
    def set_input(self, input_data):
        '''Sets the acquisition data to use for reconstruction.
        '''
        # TODO: move to C++
        assert_validity(input_data, AcquisitionData)
        _setParameter(self.handle, 'Reconstruction', \
                      'input_data', input_data.handle)
    def set_output_filename_prefix(self, prefix):
        '''Sets the output file name prefix.
        '''
        _set_char_par\
            (self.handle, 'Reconstruction', 'output_filename_prefix', prefix)
    def reconstruct(self, image):
        '''Performs reconstruction.
        '''
        assert_validity(image, ImageData)
        try_calling(pystir.cSTIR_runReconstruction(self.handle, image.handle))
    def get_output(self):
        '''Returns the reconstructed image.
        '''
        # TODO: move to C++
        return self.image

class IterativeReconstructor(Reconstructor):
    '''
    Class for a generic iterative PET reconstructor.
    '''
    def __init__(self):
        self.handle = None
        self.image = None
        self.subset = 0
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def set_num_subsets(self, n):
        '''Same as in ObjectiveFunction.
        '''
        _set_int_par\
            (self.handle, 'IterativeReconstruction', 'num_subsets', n)
    def get_num_subsets(self):
        return _int_par\
            (self.handle, 'IterativeReconstruction', 'num_subsets')
##    def set_start_subset_num(self, n):
##        _set_int_par\
##            (self.handle, 'IterativeReconstruction', 'start_subset_num', n)
##    def get_start_subset_num(self):
##        return _int_par\
##            (self.handle, 'IterativeReconstruction', 'start_subset_num')
    def set_num_subiterations(self, n):
        '''Sets number of iterations.

        Each iteration works with a subset, and is therefore referred to as
        subiteration.
        '''
        _set_int_par\
            (self.handle, 'IterativeReconstruction', 'num_subiterations', n)
    def get_num_subiterations(self):
        '''Returns the number of subiterations.
        '''
        return _int_par\
            (self.handle, 'IterativeReconstruction', 'num_subiterations')
##    def set_start_subiteration_num(self, n):
##        _set_int_par\
##            (self.handle, 'IterativeReconstruction', 'start_subiteration_num', n)
##    def get_start_subiteration_num(self):
##        return _int_par\
##            (self.handle, 'IterativeReconstruction', 'start_subiteration_num')
##    def set_subiteration_num(self, iter):
##        _set_int_par\
##            (self.handle, 'IterativeReconstruction', 'subiteration_num', iter)
    def get_subiteration_num(self):
        '''Returns the number of current subiteration.
        '''
        return _int_par\
            (self.handle, 'IterativeReconstruction', 'subiteration_num')
    def set_save_interval(self, n):
        '''Defines how often to save image iterates (n = 1: on each subiteration
        etc.)
        '''
        _set_int_par\
            (self.handle, 'IterativeReconstruction', 'save_interval', n)
##    def set_inter_iteration_filter_interval(self, n):
##        _set_int_par\
##            (self.handle, 'IterativeReconstruction',\
##             'inter_iteration_filter_interval', n)
    def set_objective_function(self, obj):
        '''Specifies the objective function to be maximized.
        '''
        assert_validity(obj, ObjectiveFunction)
        _setParameter\
            (self.handle, 'IterativeReconstruction',\
             'objective_function', obj.handle)
##    def get_objective_function(self):
##        obj_fun = ObjectiveFunction()
##        obj_fun.handle = pystir.cSTIR_parameter\
##            (self.handle, 'IterativeReconstruction', 'objective_function')
##        check_status(obj_fun.handle)
##        return obj_fun
##    def set_inter_iteration_filter(self, f):
##        pystir.cSTIR_setParameter\
##            (self.handle, 'IterativeReconstruction',\
##             'inter_iteration_filter_type', f.handle)
##    def get_inter_iteration_filter(self):
##        filter = DataProcessor()
##        filter.handle = pystir.cSTIR_parameter\
##            (self.handle, 'IterativeReconstruction',\
##             'inter_iteration_filter_type')
##        check_status(filter.handle)
##        return filter
    def set_up(self, image):
        '''Sets up the object.
        '''
        assert_validity(image, ImageData)
        handle = pystir.cSTIR_setupReconstruction(self.handle, image.handle)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
    def set_current_estimate(self, image):
        '''Sets image estimate.
        '''
        assert_validity(image, ImageData)
        self.image = image
    def process(self):
        '''Performs reconstruction.
        '''
        if self.image is None:
            raise error('current estimate not set')
        try_calling(pystir.cSTIR_runReconstruction\
                    (self.handle, self.image.handle))
    def get_current_estimate(self):
        '''Returns current image estimate.
        '''
        return self.image
    def update_current_estimate(self):
        '''Updates current image estimate by performing one subiteration.
        '''
        if self.image is None:
            raise error('current estimate not set')
        assert_validity(self.image, ImageData)
        try_calling(pystir.cSTIR_updateReconstruction\
                    (self.handle, self.image.handle))
    def set_current_subset_num(self, subset):
        '''Sets subset to be used on the next iteration..
        '''
        self.subset = subset
    def get_subset_sensitivity(self):
        '''
        Returns an ImageData object containing sensitivity image for the
        specified subset.
        '''
        obj_fun = self.get_objective_function()
        return obj_fun.get_subset_sensitivity(self.subset)
    def update(self, image):
        '''Applies one subiteration to the image estimate passed as the
        argument.
        '''
        assert_validity(image, ImageData)
        self.set_current_estimate(image)
        self.update_current_estimate()
        return self.get_current_estimate()

class OSMAPOSLReconstructor(IterativeReconstructor):
    '''
    Class for reconstructor objects using Ordered Subsets Maximum A Posteriori 
    One Step Late reconstruction algorithm, see
    http://stir.sourceforge.net/documentation/doxy/html/classstir_1_1OSMAPOSLReconstruction.html
    '''
    def __init__(self, filename = ''):
        self.handle = None
        self.image = None
        self.name = 'OSMAPOSL'
        self.handle = pystir.cSTIR_objectFromFile\
            ('OSMAPOSLReconstruction', filename)
        check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
##    def set_MAP_model(self, model):
##        _set_char_par\
##            (self.handle, self.name, 'MAP_model', model)
##    def get_objective_function(self):
##        obj_fun = PoissonLogLikelihoodWithLinearModelForMean()
##        obj_fun.handle = pystir.cSTIR_parameter\
##            (self.handle, self.name, 'objective_function')
##        check_status(obj_fun.handle)
##        return obj_fun

class OSSPSReconstructor(IterativeReconstructor):
    '''
    Class for reconstructor objects using Ordered Subsets Separable 
    Paraboloidal Surrogate reconstruction algorithm, see
    http://stir.sourceforge.net/documentation/doxy/html/classstir_1_1OSSPSReconstruction.html
    '''
    def __init__(self, filename = ''):
        self.handle = None
        self.image = None
        self.name = 'OSSPS'
        self.handle = pystir.cSTIR_objectFromFile\
                      ('OSSPSReconstruction', filename)
        check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def set_relaxation_parameter(self, value):
        _set_float_par\
            (self.handle, self.name, 'relaxation_parameter', value)

def make_Poisson_loglikelihood(acq_data, model = 'LinearModelForMean'):
    '''
    Selects the objective function based on the acquisition data and acquisition
    model types.
    '''
    # only this objective function is implemented for now
    obj_fun = PoissonLogLikelihoodWithLinearModelForMeanAndProjData()
    obj_fun.set_acquisition_data(acq_data)
    return obj_fun
