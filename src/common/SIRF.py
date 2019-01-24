''' 
Object-Oriented wrap for the cSIRF-to-Python interface pysirf.py
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
import numpy
try:
    import pylab
    HAVE_PYLAB = True
except:
    HAVE_PYLAB = False
import sys

from sirf.Utilities import assert_validities, check_status
import pyiutilities as pyiutil
import sirf.pysirf as pysirf

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
    def multiply(self, other):
        '''
        Returns the elementwise product of this and another container 
        data viewed as vectors.
        other: DataContainer
        '''
        assert_validities(self, other)
        z = self.same_object()
        z.handle = pysirf.cSIRF_multiply(self.handle, other.handle)
        check_status(z.handle)
        return z
    def divide(self, other):
        '''
        Returns the elementwise ratio of this and another container 
        data viewed as vectors.
        other: DataContainer
        '''
        assert_validities(self, other)
        z = self.same_object()
        z.handle = pysirf.cSIRF_divide(self.handle, other.handle)
        check_status(z.handle)
        return z
    def __add__(self, other):
        '''
        Overloads + for data containers.

        Returns the sum of the container data with another container 
        data viewed as vectors.
        other: DataContainer
        '''
        assert_validities(self, other)
        one = numpy.asarray([1.0, 0.0], dtype = numpy.float32)
        z = self.same_object()
        z.handle = pysirf.cSIRF_axpby \
            (one.ctypes.data, self.handle, one.ctypes.data, other.handle)
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
        pl_one = numpy.asarray([1.0, 0.0], dtype = numpy.float32)
        mn_one = numpy.asarray([-1.0, 0.0], dtype = numpy.float32)
        z = self.same_object()
        z.handle = pysirf.cSIRF_axpby \
            (pl_one.ctypes.data, self.handle, mn_one.ctypes.data, other.handle)
        check_status(z.handle)
        return z;
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
    def __truediv__(self, other):
        '''
        Overloads / for data containers multiplication by a scalar or another
        data container.

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

class ImageData(DataContainer):
    def write(self, filename):
        """Write image to file."""
    pass
