''' 
Object-Oriented wrap for the cSIRF-to-Python interface pysirf.py
'''

## SyneRBI Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2020 Rutherford Appleton Laboratory STFC
## Copyright 2015 - 2020 University College London
##
## This is software developed for the Collaborative Computational
## Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
## (http://www.ccpsynerbi.ac.uk/).
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
import warnings

from sirf.Utilities import assert_validity, assert_validities, \
     cpp_int_dtype, check_status, try_calling, error
import pyiutilities as pyiutil
import sirf.pysirf as pysirf

from numbers import Number
import deprecation

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
        print("SIRF.DataContainer __del__ with handle {}.".format(self.handle))
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def __neg__(self):
        return self.multiply(-1.0)

    def __add__(self, other):
        '''
        Overloads + for data containers.

        Returns the sum of the container data with another container 
        data viewed as vectors.
        other: DataContainer
        '''
        return self.add(other)

    def __sub__(self, other):
        '''
        Overloads - for data containers.

        Returns the difference of the container data with another container 
        data viewed as vectors.
        other: DataContainer
        '''
        if isinstance(other, (DataContainer, Number) ):
            return self.subtract(other)
        return NotImplemented

    def __mul__(self, other):
        '''
        Overloads * for data containers multiplication by a scalar or another
        data container.

        Returns the product self*other if other is a scalar
        or the elementwise product if other is of the same type as self.
        other: DataContainer or a (real or complex) scalar
        '''
        return self.multiply(other)

    def __rmul__(self, other):
        '''
        Overloads * for data containers multiplication by a scalar from
        the left, i.e. computes and returns the product other*self.
        other: a real or complex scalar
        '''
        return self * other

    def __truediv__(self, other):
        '''
        Overloads / for data containers division by a scalar or (elementwise)
        another data container (Python 3.*)

        Returns the ratio self/other if other is a scalar
        or the elementwise ratio if other is of the same type as self.
        other: DataContainer or a (real or complex) scalar
        '''
        return self.divide(other)

    def __iadd__(self, other):
        self.add(other, out=self)
        return self

    def __imul__(self, other):
        self.multiply(other, out=self)
        return self

    def __isub__(self, other):
        self.subtract(other, out=self)
        return self

    def __itruediv__(self, other):
        self.divide(other, out=self)
        return self

    def same_object(self):
        '''
        Returns an object of the same type as self.

        Since this class is abstract, its methods cannot itself create a new
        object when e.g. adding two objects of this class, so new object is
        created by the first object using its same_object() method - see
        clone() below.
        '''
        return DataContainer()

    def clone(self):
        assert self.handle is not None
        x = self.same_object()
        x.handle = pysirf.cSIRF_clone(self.handle)
        check_status(x.handle)
        return x

    def copy(self):
        '''alias of clone'''
        return self.clone()

    def get_uniform_copy(self, value=1.0):
        '''Initialises an instance of DataContainer based on the template'''
        y = self.clone()
        y.fill(value)
        return y

    def allocate(self, value=0, **kwargs):
        """Allocates a copy of self and fills with values

        value: Python float or str
            float: the value to fill with
            'random': fill with random values ranging between 0 1nd 1 generated
                by numpy.random.random_sample, optionally using seed provided by
                kwarg 'seed'
            'random_int': fill with random integers ranging between 0 and the
                value optionally provided by kwarg 'max_value' (by default, 100)

        CIL/SIRF compatibility
        """
        if value in ['random', 'random_int']:
            out = self.get_uniform_copy()
            shape = out.shape
            seed = kwargs.get('seed', None)
            if seed is not None:
                numpy.random.seed(seed)
            if value == 'random':
                out.fill(numpy.random.random_sample(shape))
            elif value == 'random_int':
                max_value = kwargs.get('max_value', 100)
                out.fill(numpy.random.randint(max_value,size=shape))
        else:
            out = self.get_uniform_copy(value)
        return out

    def write(self, filename):
        '''
        Writes to file.
        '''
        try_calling(pysirf.cSIRF_write(self.handle, filename))

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

    def is_empty(self):
        return self.number() < 1

    def is_complex(self):
        assert self.handle is not None
        handle = pysirf.cSIRF_isComplex(self.handle)
        check_status(handle)
        i = pyiutil.intDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return i != 0

    def conjugate(self, out=None):
        ''' Computes complex conjugate of self.

            Use y = x.conjugate() to get the conjugated copy of x.
            Use x.conjugate(out=x) to conjugate in-place.
        '''
        if out is self:
            try_calling(pysirf.cSIRF_conjugate(self.handle))
            return
        elif out is None:
            x = self.same_object()
        else:
            x = out
        if x.handle is not None:
            pyiutil.deleteDataHandle(x.handle)
        x.handle = pysirf.cSIRF_conjugated(self.handle)
        check_status(x.handle)
        if out is None:
            return x

    def norm(self):
        '''
        Returns the 2-norm of the container data viewed as a vector.
        '''
        handle = pysirf.cSIRF_norm(self.handle)
        check_status(handle)
        r = pyiutil.floatDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return r;

    def squared_norm(self):
        '''Returns the squared norm of a DataContainer viewed as a vector

        CIL/SIRF compatibility'''
        return self.norm() ** 2

    def dot(self, other):
        '''
        Returns the dot product of the container data with another container 
        data viewed as vectors.
        other: DataContainer
        '''
        assert_validities(self, other)
        # Check if input are the same size
        if self.size != other.size:
            raise ValueError("Input sizes are expected to be equal, got " \
                + self.size + " and " + other.size + " instead.")
        z = numpy.zeros((2,), dtype=numpy.float32)
        try_calling(pysirf.cSIRF_compute_dot(self.handle, other.handle, z.ctypes.data))
        if z[1] == 0:
            return z[0]
        return z[0] + 1j*z[1]

    def sum(self):
        '''
        Returns the sum of the elements of self data
        '''
        z = numpy.zeros((2,), dtype=numpy.float32)
        try_calling(pysirf.cSIRF_compute_sum(self.handle, z.ctypes.data))
        if z[1] == 0:
            return z[0]
        return z[0] + 1j*z[1]

    def max(self):
        '''
        Returns the maximum of the elements of self data
        '''
        z = numpy.zeros((2,), dtype=numpy.float32)
        try_calling(pysirf.cSIRF_compute_max(self.handle, z.ctypes.data))
        if z[1] == 0:
            return z[0]
        return z[0] + 1j*z[1]

    def add(self, other, out=None):
        '''
        Addition for data containers.

        If other is a DataContainer, returns the sum of data
        stored in self and other viewed as vectors.
        If other is a scalar, returns the same with the second vector filled
        with the value of other.
        other: DataContainer or scalar.
        out:   DataContainer to store the result to.
        '''
        if out is None:
            z = self.same_object()
        else:
            z = out
            assert_validities(self, z)
        if isinstance(other, Number):
            a = numpy.asarray([other.real, other.imag], dtype=numpy.float32)
            if out is None:
                z.handle = pysirf.cSIRF_sum(self.handle, a.ctypes.data)
                check_status(z.handle)
            else:
                try_calling(pysirf.cSIRF_add(self.handle, a.ctypes.data, z.handle))
        else:
            assert_validities(self, other)
            one = numpy.asarray([1.0, 0.0], dtype = numpy.float32)
            if out is None:
                z.handle = pysirf.cSIRF_axpby(one.ctypes.data, self.handle, one.ctypes.data, other.handle)
                check_status(z.handle)
            else:
                try_calling(pysirf.cSIRF_axpbyAlt(one.ctypes.data, self.handle, one.ctypes.data, other.handle, z.handle))
        return z

    def subtract(self, other, out=None):
        '''
        Subtraction for data containers.

        If other is a DataContainer, returns the difference of data
        stored in self and other viewed as vectors.
        If other is a scalar, returns the same with the second vector filled
        with the value of other.
        other: DataContainer or scalar.
        other: DataContainer
        '''
        if not isinstance (other, (DataContainer, Number)):
            return NotImplemented
        if isinstance(other, Number):
            return self.add(-other, out=out)
        assert_validities(self, other)
        pl_one = numpy.asarray([1.0, 0.0], dtype = numpy.float32)
        mn_one = numpy.asarray([-1.0, 0.0], dtype = numpy.float32)
        if out is None:
            z = self.same_object()
            z.handle = pysirf.cSIRF_axpby \
                (pl_one.ctypes.data, self.handle, mn_one.ctypes.data, other.handle)
            check_status(z.handle)
            return z
        else:
            assert_validities(self, out)
            z = out
            try_calling(pysirf.cSIRF_axpbyAlt \
                (pl_one.ctypes.data, self.handle, mn_one.ctypes.data, other.handle, z.handle))

    def multiply(self, other, out=None):
        '''
        Elementwise multiplication for data containers.

        Returns the elementwise product of data stored in self and other.
        other: DataContainer or scalar.
        out:   DataContainer to store the result to.
        '''
        if out is None:
            return self.binary(other, 'multiply')
        self.binary(other, 'multiply', out=out)

    def divide(self, other, out=None):
        '''
        Elementwise ratio for data containers.

        Returns the elementwise ratio of data stored in self and other.
        other: DataContainer or scalar.
        out:   DataContainer to store the result to.
        '''
        if isinstance(other, Number):
            if out is None:
                return self.binary(1./other, 'multiply')
            self.binary(1./other, 'multiply', out=out)
        else:
            if out is None:
                return self.binary(other, 'divide')
            self.binary(other, 'divide', out=out)

    def maximum(self, other, out=None):
        '''
        Elementwise maximum for data containers.

        Returns the elementwise maximum of data stored in self and other.
        other: DataContainer or scalar.
        out:   DataContainer to store the result to.
        '''
        if out is None:
            return self.binary(other, 'maximum')
        self.binary(other, 'maximum', out=out)

    def minimum(self, other, out=None):
        '''
        Elementwise minimum for data containers.

        Returns the elementwise minimum of data stored in self and other.
        other: DataContainer or scalar.
        out:   DataContainer to store the result to.
        '''
        if out is None:
            return self.binary(other, 'minimum')
        self.binary(other, 'minimum', out=out)

    @deprecation.deprecated(details="Please use the sapyb method instead")
    def axpby(self, a, b, y, out=None, **kwargs):
        '''
        Linear combination for data containers.

        Returns the linear combination of the self data with another container 
        data y viewed as vectors.
        a: multiplier to self, can be a number or a DataContainer
        b: multiplier to y, can be a number or a DataContainer 
        y: DataContainer
        out:   DataContainer to store the result to.
        '''
        return self.sapyb(a, y, b, out=out, **kwargs)

    def sapyb(self, a, y, b, out=None, **kwargs):
        '''
        Linear combination for data containers: new interface.

        Returns the linear combination of the self data with another container 
        data y viewed as vectors.
        a: multiplier to self, can be a number or a DataContainer
        b: multiplier to y, can be a number or a DataContainer 
        y: DataContainer
        out:   DataContainer to store the result to, can be self or y.
        '''

        assert_validities(self, y)

        if out is not None:
            assert_validities(self, out)
            z = out
        else:
            z = self.same_object()

        if isinstance(a, Number):
            alpha = numpy.asarray([a.real, a.imag], dtype = numpy.float32)
            if isinstance(b, Number):
                #a is scalar, b is scalar
                beta = numpy.asarray([b.real, b.imag], dtype = numpy.float32)
                if out is None:
                    z.handle = pysirf.cSIRF_axpby(alpha.ctypes.data, self.handle, beta.ctypes.data, y.handle)
                else:
                    try_calling(pysirf.cSIRF_axpbyAlt(alpha.ctypes.data, self.handle, beta.ctypes.data, y.handle, z.handle))
            else:
                #a is scalar, b is array
                if out is None:
                    z.handle = pysirf.cSIRF_XapYB(self.handle, alpha.ctypes.data, y.handle, b.handle)
                else:
                    try_calling(pysirf.cSIRF_XapYBAlt(self.handle, alpha.ctypes.data, y.handle, b.handle, z.handle))
        else:
            assert_validities(self, a)
            if isinstance(b, Number):
                #a is array, b is scalar
                beta = numpy.asarray([b.real, b.imag], dtype = numpy.float32)
                if out is None:
                    z.handle = pysirf.cSIRF_XapYB(y.handle, beta.ctypes.data, self.handle, a.handle)
                else:
                    try_calling(pysirf.cSIRF_XapYBAlt(y.handle, beta.ctypes.data, self.handle, a.handle, z.handle))
            else:
                #a is array, b is array
                assert_validities(self, b)
                if out is None:
                    z.handle = pysirf.cSIRF_xapyb(self.handle, a.handle, y.handle, b.handle)
                else:
                    try_calling(pysirf.cSIRF_xapybAlt(self.handle, a.handle, y.handle, b.handle, z.handle))

        if out is None:
            check_status(z.handle)
            return z

    def power(self, other, out=None):
        '''Power function for DataContainers

        SIRF/CIL compatibility
        '''
        if out is None:
            return self.binary(other, 'power')
        self.binary(other, 'power', out=out)

    def abs(self, out=None):
        '''Returns the element-wise absolute value of the DataContainer data

        '''
        if out is None:
            return self.unary('abs')
        self.unary('abs', out=out)

    def sign(self, out=None):
        '''Returns the element-wise sign of the DataContainer data

        '''
        if out is None:
            return self.unary('sign')
        self.unary('sign', out=out)

    def sqrt(self, out=None):
        '''Returns the element-wise sqrt of the DataContainer data

        '''
        if out is None:
            return self.unary('sqrt')
        self.unary('sqrt', out=out)

    def exp(self, out=None):
        '''Returns the element-wise exp of the DataContainer data

        '''
        if out is None:
            return self.unary('exp')
        self.unary('exp', out=out)

    def log(self, out=None):
        '''Returns the element-wise log of the DataContainer data

        '''
        if out is None:
            return self.unary('log')
        self.unary('log', out=out)

    def binary(self, other, f, out=None):
        '''Applies function f(x,y) element-wise to self and other.

        other: DataContainer or Number
        f: the name of the function to apply, Python str.
        '''
        if out is None:
            out = self.same_object()
        if isinstance(other, Number):
            y = numpy.asarray([other.real, other.imag], dtype=numpy.float32)
            if out.handle is None:
                out.handle = pysirf.cSIRF_semibinary(self.handle, y.ctypes.data, f)
                check_status(out.handle)
                return out
            else:
                try_calling(pysirf.cSIRF_compute_semibinary(self.handle, y.ctypes.data, \
                                                        f, out.handle))
        else:
            assert_validities(self, other)
            if out.handle is None:
                out.handle = pysirf.cSIRF_binary(self.handle, other.handle, f)
                check_status(out.handle)
                return out
            else:
                try_calling(pysirf.cSIRF_compute_binary(self.handle, other.handle, \
                                                        f, out.handle))

    def unary(self, f, out=None):
        '''Applies function f(x) element-wise to self data.

        f: the name of the function to apply, Python str.
        '''
        if out is None:
            out = self.same_object()
        if out.handle is None:
            out.handle = pysirf.cSIRF_unary(self.handle, f)
            check_status(out.handle)
            return out
        else:
            assert_validities(self, out)
            try_calling(pysirf.cSIRF_compute_unary(self.handle, f, out.handle))

    @property
    def shape(self):
        '''Returns the shape of the object data
        '''
        if self.is_empty():
            return (0,)
        return self.dimensions()

    @property
    def size(self):
        '''Returns the number of elements in the object data.'''
        if self.is_empty():
            return 0
        return numpy.prod(self.dimensions())

    @property
    def dtype(self):
        handle = pysirf.cSIRF_bits(self.handle)
        check_status(handle)
        bits = pyiutil.intDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        if self.is_complex():
            dt = 'complex%s' % bits
        else:
            dt = 'float%s' % bits
        return numpy.dtype(dt)
    

class ImageData(DataContainer):
    '''
    Image data ABC
    '''

    def equal(self, other):
        '''
        Overloads == for ImageData.

        other: ImageData
        '''
        assert_validity(self, ImageData)
        assert_validity(other, ImageData)
        handle = pysirf.cSIRF_equalImages(self.handle, other.handle)
        check_status(handle)
        same = pyiutil.intDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return same

    def __eq__(self, other):
        return self.equal(other)

    def __ne__(self, other):
        '''
        Overloads != for ImageData.

        other: ImageData
        '''
        return not (self == other)

    def read(self, file, engine, verb):
        self.handle = pysirf.cSIRF_readImageData(file, engine, verb)
        check_status(self.handle)

    def fill(self, image):
        try_calling(pysirf.cSIRF_fillImageFromImage(self.handle, image.handle))
        return self

    def get_geometrical_info(self):
        """Get the image's geometrical info."""
        try:
            import sirf.STIR
            if isinstance(self, sirf.STIR.ImageData):
                warnings.warn("geometrical info for STIR.ImageData might be incorrect")
        except:
            pass
        geom_info = GeometricalInfo()
        geom_info.handle = pysirf.cSIRF_ImageData_get_geom_info(self.handle)
        check_status(geom_info.handle)
        return geom_info

    def reorient(self, geom_info):
        """Reorient image. Requires that dimensions match."""
        if not isinstance(geom_info, GeometricalInfo):
            raise AssertionError()
        try_calling(pysirf.cSIRF_ImageData_reorient(self.handle, geom_info.handle))

DataContainer.register(ImageData)

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

    @deprecation.deprecated(details="Please use get_info method instead")
    def print_info(self):
        """Print the geom info"""
        print(self.get_info())

    def get_info(self):
        """Return the geom info as string"""
        handle = pysirf.cSIRF_GeomInfo_get(self.handle)
        check_status(handle)
        info = pyiutil.charDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return info

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
        arr = numpy.ndarray((3,), dtype = cpp_int_dtype())
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
