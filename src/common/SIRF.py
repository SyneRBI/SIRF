''' 
Object-Oriented wrap for the cSIRF-to-Python interface pysirf.py
'''

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
## Copyright 2015 - 2019 University College London
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
try:
    import pylab
    HAVE_PYLAB = True
except:
    HAVE_PYLAB = False
import sys

from sirf.Utilities import assert_validities, check_status, try_calling
import pyiutilities as pyiutil
import sirf.pysirf as pysirf

from numbers import Number

try:
    input = raw_input
except NameError:
    pass

if sys.version_info[0] >= 3 and sys.version_info[1] >= 4:
    ABC = abc.ABC
else:
    ABC = abc.ABCMeta('ABC', (), {})

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
    def clone(self):
        assert self.handle is not None
        x = self.same_object()
        x.handle = pysirf.cSIRF_clone(self.handle)
        check_status(x.handle)
        return x
    def number(self):
        '''
        Returns the number of items in the container.
        '''
        assert self.handle is not None
        handle = pysirf.cSIRF_dataItems(self.handle)
        check_status(handle)
        n = pyiutil.intDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return n
    def norm(self):
        '''
        Returns the 2-norm of the container data viewed as a vector.
        '''
        assert self.handle is not None
        handle = pysirf.cSIRF_norm(self.handle)
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
        handle = pysirf.cSIRF_dot(self.handle, other.handle)
        check_status(handle)
        r = pyiutil.floatDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return r
    def multiply(self, other, out=None):
        '''
        Returns the elementwise product of this and another container 
        data viewed as vectors.
        other: DataContainer
        out:   DataContainer to store the result to.
        '''
        if isinstance(other , ( Number, int, float, numpy.float32 )):
            tmp = other + numpy.zeros(self.as_array().shape)
            other = self.copy()
            other.fill(tmp)
        assert_validities(self, other)
        if out is None:
            z = self.same_object()
        else:
            assert_validities(self, out)
            z = out
        z.handle = pysirf.cSIRF_multiply(self.handle, other.handle)
        check_status(z.handle)
        return z
    def divide(self, other, out=None):
        '''
        Returns the elementwise ratio of this and another container 
        data viewed as vectors.
        other: DataContainer
        out:   DataContainer to store the result to.
        '''
        if isinstance(other , ( Number, int, float, numpy.float32 )):
            tmp = other + numpy.zeros_like(self.as_array().shape)
            other = self.copy()
            other.fill(tmp)
        assert_validities(self, other)
        if out is None:
            z = self.same_object()
        else:
            assert_validities(self, out)
            z = out
        z.handle = pysirf.cSIRF_divide(self.handle, other.handle)
        check_status(z.handle)
        return z
    def add(self, other, out=None):
        '''
        Addition for data containers.

        Returns the sum of the container data with another container 
        data viewed as vectors.
        other: DataContainer
        out:   DataContainer to store the result to.
        '''
        if isinstance(other , ( Number, int, float, numpy.float32 )):
            tmp = other + numpy.zeros(self.as_array().shape)
            other = self.copy()
            other.fill(tmp)
        assert_validities(self, other)
        one = numpy.asarray([1.0, 0.0], dtype = numpy.float32)
        if out is None:
            z = self.same_object()
        else:
            assert_validities(self, out)
            z = out
        z.handle = pysirf.cSIRF_axpby \
            (one.ctypes.data, self.handle, one.ctypes.data, other.handle)
        check_status(z.handle)
        return z;
    def write(self, filename):
        '''
        Writes to file.
        '''
        assert self.handle is not None
        try_calling(pysirf.cSIRF_write(self.handle, filename))
    def __add__(self, other):
        '''
        Overloads + for data containers.

        Returns the sum of the container data with another container 
        data viewed as vectors.
        other: DataContainer
        '''
        return self.add(other)
    def subtract(self, other, out=None):
        '''
        Overloads - for data containers.

        Returns the difference of the container data with another container 
        data viewed as vectors.
        other: DataContainer
        '''
        if isinstance(other , ( Number, int, float, numpy.float32 )):
            tmp = other + numpy.zeros(self.as_array().shape)
            other = self.copy()
            other.fill(tmp)
        assert_validities(self, other)
        pl_one = numpy.asarray([1.0, 0.0], dtype = numpy.float32)
        mn_one = numpy.asarray([-1.0, 0.0], dtype = numpy.float32)
        if out is None:
            z = self.same_object()
        else:
            assert_validities(self, out)
            z = out
        z.handle = pysirf.cSIRF_axpby \
            (pl_one.ctypes.data, self.handle, mn_one.ctypes.data, other.handle)
        check_status(z.handle)
        return z;
    def __sub__(self, other):
        '''
        Overloads - for data containers.

        Returns the difference of the container data with another container 
        data viewed as vectors.
        other: DataContainer
        '''
        return self.subtract(other)
    def __mul__(self, other):
        '''
        Overloads * for data containers multiplication by a scalar or another
        data container.

        Returns the product self*other if other is a scalar
        or the elementwise product if it is DataContainer.
        other: DataContainer or a (real or complex) scalar
        '''
        assert self.handle is not None
        if type(self) == type(other):
            return self.multiply(other)
        z = self.same_object()
        try:
            a = numpy.asarray([other.real, other.imag], dtype = numpy.float32)
            zero = numpy.zeros((2,), dtype = numpy.float32)
            z.handle = pysirf.cSIRF_axpby \
                (a.ctypes.data, self.handle, zero.ctypes.data, self.handle)
            z.src = 'mult'
            check_status(z.handle)
            return z;
        except:
            raise error('wrong multiplier')

    def __rmul__(self, other):
        '''
        Overloads * for data containers multiplication by a scalar from
        the left, i.e. computes and returns the product other*self.
        other: a real or complex scalar
        '''
        assert self.handle is not None
        z = self.same_object()
        try:
            a = numpy.asarray([other.real, other.imag], dtype = numpy.float32)
            zero = numpy.zeros((2,), dtype = numpy.float32)
            z.handle = pysirf.cSIRF_axpby \
                (a.ctypes.data, self.handle, zero.ctypes.data, self.handle)
            check_status(z.handle)
            return z;
        except:
            raise error('wrong multiplier')
    def __div__(self, other):
        '''
        Overloads / for data containers division by a scalar or (elementwise)
        another data container (Python 2.*)

        Returns the product self*other if other is a scalar
        or the elementwise product if it is DataContainer.
        other: DataContainer or a (real or complex) scalar
        '''
        assert self.handle is not None
        if type(self) == type(other):
            return self.divide(other)
        z = self.same_object()
        try:
            other = 1.0/other
            a = numpy.asarray([other.real, other.imag], dtype = numpy.float32)
            zero = numpy.zeros((2,), dtype = numpy.float32)
            z.handle = pysirf.cSIRF_axpby \
                (a.ctypes.data, self.handle, zero.ctypes.data, self.handle)
            check_status(z.handle)
            return z;
        except:
            raise error('wrong multiplier')
    def copy(self):
        '''alias of clone'''
        return self.clone()
    def power(self, other, out=None):
        '''Power function for DataContainers

        uses NumPy
        SIRF/CIL compatibility
        '''
        if out is None:
            z = self.clone()
        else:
            assert_validities(self, out)
            z = out
        
        if isinstance(other, Number):
            tmp =  numpy.power(self.as_array(), other)
            z.fill( tmp )
        else:
            assert_validities(self, other)
            z.fill(
               numpy.power(self.as_array(), other.as_array())
            )
        return z
    def maximum(self, other, out=None):
        '''Element-wise maximum of DataContainer elements.

        Compare two DataContainers and returns a new array containing the element-wise maxima. Output can be pre-allocated in variable out.

        uses NumPy
        SIRF/CIL compatibility
        '''
        if out is None:
            z = self.clone()
        else:
            assert_validities(self, out)
            z = out
        if isinstance(other, Number):
            z.fill(
               numpy.maximum(self.as_array(), other)
            )
        else:
            assert_validities(self, other)
            z.fill(
               numpy.maximum(self.as_array(), other.as_array())
            )
        return z
    def minimum(self, other, out=None):
        '''Element-wise minimum of DataContainer elements.

        Compare two DataContainers and returns a new array containing the element-wise minima. Output can be pre-allocated in variable out.

        uses NumPy
        SIRF/CIL compatibility
        '''

        if out is None:
            z = self.clone()
        else:
            assert_validities(self, out)
            z = out
        if isinstance(other, Number):
            z.fill(
               numpy.minimum(self.as_array(), other)
            )
        else:
            assert_validities(self, other)
            z.fill(
               numpy.minimum(self.as_array(), other.as_array())
            )
        return z
    # inline algebra
    def __iadd__(self, other):
        '''Not quite in-place add'''
        self.fill(self.add(other))
        return self
    def __imul__(self, other):
        '''Not quite in-place multiplication'''
        if isinstance(other, Number):
            z = other * self
            self.fill(z.as_array())
            return self
        self.fill(self.multiply(other).as_array())
        return self
    def __isub__(self, other):
        '''Not quite in-place subtract'''
        self.fill(self.subtract(other).as_array())
        return self
    def __idiv__(self, other):
        '''Not quite in-place division'''
        if isinstance(other, Number):
            z = (1./other) * self
            self.fill(z.as_array())
            return self
        self.fill(self.divide(other).as_array())
        return self
    def abs(self, out=None):
        '''Returns the element-wise absolute value of the DataContainer data
        
           uses NumPy 
        '''
        if out is None:
            z = self.clone()
        else:
            assert_validities(self, out)
            z = out
        z.fill(
               numpy.abs(self.as_array())
        )
        return z
    def sign(self, out=None):
        '''Returns the element-wise sign of the DataContainer data
        
           uses NumPy 
        '''
        if out is None:
            z = self.clone()
        else:
            assert_validities(self, out)
            z = out
        z.fill(
               numpy.sign(self.as_array())
        )
        return z
    def sqrt(self, out=None):
        '''Returns the element-wise sqrt of the DataContainer data

           uses NumPy
        '''
        if out is None:
            z = self.clone()
        else:
            assert_validities(self, out)
            z = out
        z.fill(
               numpy.sqrt(self.as_array())
        )
        return z
    def exp(self, out=None):
        '''Returns the element-wise exp of the DataContainer data

           uses NumPy
        '''
        if out is None:
            z = self.clone()
        else:
            assert_validities(self, out)
            z = out
        z.fill(
               numpy.exp(self.as_array())
        )
        return z
    def log(self, out=None):
        '''Returns the element-wise log of the DataContainer data

           uses NumPy
        '''
        if out is None:
            z = self.clone()
        else:
            assert_validities(self, out)
            z = out
        z.fill(
               numpy.log(self.as_array())
        )
        return z
    def sum(self):
        '''Returns the sum of DataContainer elements.

           it is a reduction operation

           uses NumPy
        '''
        return numpy.sum(self.as_array())
    def get_uniform_copy(self, value=1.0):
        '''Initialises an instance of DataContainer based on the template'''
        y = self.clone()
        y.fill(value)
        return y

    def squared_norm(self):
        '''Returns the squared norm of a DataContainer viewed as a vector
        
        CIL/SIRF compatibility'''
        return self.norm() ** 2

    def __truediv__(self, other):
        '''
        Same as __div__ but for Python 3.*
        '''
        return self.__div__(other)
    @property
    def shape(self):
        '''returns the shape of the data array
        
        CIL/SIRF compatibility
        '''
        return self.as_array().shape

class ImageData(DataContainer):
    def fill(self, image):
        try_calling(pysirf.cSIRF_fillImageFromImage(self.handle, image.handle))

    def get_geometrical_info(self):
        """Get the image's geometrical info."""
        geom_info = GeometricalInfo()
        geom_info.handle = pysirf.cSIRF_ImageData_get_geom_info(self.handle)
        check_status(geom_info.handle)
        return geom_info


class DataHandleVector(object):
    """
    DataHandle vector.
    """
    def __init__(self):
        self.name = 'DataHandleVector'
        self.handle = pysirf.cSIRF_newObject(self.name)
        check_status(self.handle)

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def push_back(self, handle):
        """Push back new data handle."""
        try_calling(pysirf.cSIRF_DataHandleVector_push_back(self.handle, handle))
        check_status(self.handle)

class GeometricalInfo(object):
    """
    Get the geometrical information in LPS space. These are encoded
    as size (number of voxels), spacing, offset (distance to first voxel)
    and a direction matrix.
    """
    def __init__(self):
        self.name = 'GeometricalInfo'
        self.handle = None

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def print_info(self):
        """Print the geom info"""
        try_calling(pysirf.cSIRF_GeomInfo_print(self.handle))

    def get_offset(self):
        """Offset is the LPS coordinate of the centre of the first voxel."""
        arr = numpy.ndarray((3,), dtype = numpy.float32)
        try_calling(pysirf.cSIRF_GeomInfo_get_offset(self.handle, arr.ctypes.data))
        return tuple(arr)

    def get_spacing(self):
        """Spacing is the physical distance between voxels in each dimension."""
        arr = numpy.ndarray((3,), dtype = numpy.float32)
        try_calling (pysirf.cSIRF_GeomInfo_get_spacing(self.handle, arr.ctypes.data))
        return tuple(arr)
    
    def get_size(self):
        """Size is the number of voxels in each dimension."""
        arr = numpy.ndarray((3,), dtype = numpy.int32)
        try_calling (pysirf.cSIRF_GeomInfo_get_size(self.handle, arr.ctypes.data))
        return tuple(arr)

    def get_direction_matrix(self):
        """Each row gives a vector dictating the direction of the axis in LPS physical space."""
        arr = numpy.ndarray((3,3), dtype = numpy.float32)
        try_calling (pysirf.cSIRF_GeomInfo_get_direction_matrix(self.handle, arr.ctypes.data))
        return arr

    def get_index_to_physical_point_matrix(self):
        """Get the 4x4 affine matrix that converts an index to a point in LPS physical space."""
        arr = numpy.ndarray((4,4), dtype = numpy.float32)
        try_calling (pysirf.cSIRF_GeomInfo_get_index_to_physical_point_matrix(self.handle, arr.ctypes.data))
        return arr
