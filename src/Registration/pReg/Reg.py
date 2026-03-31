"""Object-Oriented wrap for the cReg-to-Python interface pyreg.py."""

# SyneRBI Synergistic Image Reconstruction Framework (SIRF)
# Copyright 2018 - 2021 University College London
# Copyright 2018 - 2021 Science Technology Facilities Council
#
# This is software developed for the Collaborative Computational
# Project in Synergistic Reconstruction for Biomedical Imaging
# (formerly CCP PETMR)
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
import sys
from numbers import Number

import numpy
import pyreg

import sirf.Reg_params as parms
from sirf import SIRF
from sirf.config import SIRF_HAS_SPM
from sirf.SIRF import ContiguousError
from sirf.Utilities import Handle, cpp_int_dtype, error, format_numpy_array_for_setter

if sys.version_info[0] >= 3 and sys.version_info[1] >= 4:
    ABC = abc.ABC
else:
    ABC = abc.ABCMeta('ABC', (), {})

INFO_CHANNEL = 0
WARNING_CHANNEL = 1
ERROR_CHANNEL = 2
ALL_CHANNELS = -1


class MessageRedirector:
    """Class for registration printing redirection to files/stdout/stderr."""
    def __init__(self, info=None, warn='stdout', errr='stdout'):
        """Create MessageRedirector object.

        This redirects Reg's ouput
        produced by info(), warning() and error(0 functions to destinations
        specified respectively by info, warn and err arguments.
        The argument values other than None, stdout, stderr, cout and cerr
        are interpreted as filenames.
        None and empty string value suppresses printing.
        """
        if info is None:
            info = ''
        if not isinstance(info, str):
            raise error('wrong info argument for MessageRedirector constructor')
        elif info in {'stdout', 'stderr', 'cout', 'cerr'}:
            self.info = pyreg.newTextPrinter(info)
            self.info_case = 0
        else:
            self.info = pyreg.newTextWriter(info)
            self.info_case = 1
        pyreg.openChannel(0, self.info)

        if warn is None:
            warn = ''
        if not isinstance(warn, str):
            raise error('wrong warn argument for MessageRedirector constructor')
        elif warn in {'stdout', 'stderr', 'cout', 'cerr'}:
            self.warn = pyreg.newTextPrinter(warn)
            self.warn_case = 0
        else:
            self.warn = pyreg.newTextWriter(warn)
            self.warn_case = 1
        pyreg.openChannel(1, self.warn)

        if errr is None:
            errr = ''
        if not isinstance(errr, str):
            raise error('wrong errr argument for MessageRedirector constructor')
        elif errr in {'stdout', 'stderr', 'cout', 'cerr'}:
            self.errr = pyreg.newTextPrinter(errr)
            self.errr_case = 0
        else:
            self.errr = pyreg.newTextWriter(errr)
            self.errr_case = 1
        pyreg.openChannel(2, self.errr)

    def __del__(self):
        """del."""
        if self.info_case == 0:
            Handle(pyreg.deleteTextPrinter(self.info))
        else:
            Handle(pyreg.deleteTextWriter(self.info))
        pyreg.closeChannel(0, self.info)
        if self.warn_case == 0:
            Handle(pyreg.deleteTextPrinter(self.warn))
        else:
            Handle(pyreg.deleteTextWriter(self.warn))
        pyreg.closeChannel(1, self.warn)
        if self.errr_case == 0:
            Handle(pyreg.deleteTextPrinter(self.errr))
        else:
            Handle(pyreg.deleteTextWriter(self.errr))
        pyreg.closeChannel(2, self.errr)


###########################################################


class _Transformation(ABC):
    """Abstract base class for transformations."""
    def __init__(self):
        """init."""
        self.handle = None
        self.name = 'Transformation'

    def get_as_deformation_field(self, ref):
        """Get any type of transformation as a deformation field.

        This is useful for joining them together. Require a reference
        image for converting transformation matrices to deformations.
        """
        if not isinstance(ref, NiftiImageData3D):
            raise AssertionError()
        output = NiftiImageData3DDeformation()
        output.handle = self.handle.cReg_Transformation_get_as_deformation_field(self.name, ref)
        return output


class NiftiImageData(SIRF.ImageData):
    """General class for nifti image data."""
    _ADD = 0
    _SUBTRACT = 1
    _MULTIPLY = 2
    _DIVIDE = 3

    def __init__(self, src=None):
        """init."""
        self.handle = None
        self.name = 'NiftiImageData'
        if src is None:
            self.handle = Handle(None, -1).cReg_newObject(self.name)
        elif isinstance(src, str):
            self.handle = Handle(None, -1).cReg_objectFromFile(self.name, src)
        elif isinstance(src, SIRF.ImageData):
            # src is ImageData
            dim = src.dimensions()
            self.handle = src.handle.cReg_NiftiImageData_from_SIRFImageData(int(dim[0] > 1))
        else:
            raise error('Wrong source in NiftiImageData constructor')

    def equal(self, other):
        """Overload comparison operator."""
        if not isinstance(other, NiftiImageData):
            raise AssertionError()
        if not self.handle.valid or not other.handle.valid:
            error('Cannot compare images as at least one is uninitialised')
        h = self.handle.cReg_NiftiImageData_equal(other)
        return int(h)

    def __eq__(self, other):
        """Overload == operator."""
        try:
            return self.equal(other)
        except:
            return super().equal(other)

    def __ne__(self, other):
        """Overload comparison operator."""
        return not self == other

    def write(self, filename, datatype=-1):
        """Save to file.

        See nifti1.h for datatypes (e.g., float (NIFTI_TYPE_FLOAT32) = 16).
        Image's original datatpye is used by default.
        """
        self.handle.cReg_NiftiImageData_write(filename, datatype)

    def get_max(self):
        """Get max."""
        return parms.float_par(self.handle, 'NiftiImageData', 'max')

    def get_min(self):
        """Get min."""
        return parms.float_par(self.handle, 'NiftiImageData', 'min')

    def get_mean(self):
        """Get mean."""
        return parms.float_par(self.handle, 'NiftiImageData', 'mean')

    def get_variance(self):
        """Get variance."""
        return parms.float_par(self.handle, 'NiftiImageData', 'variance')

    def get_standard_deviation(self):
        """Get standard deviation."""
        return parms.float_par(self.handle, 'NiftiImageData', 'std')

    def get_sum(self):
        """Get sum."""
        return parms.float_par(self.handle, 'NiftiImageData', 'sum')

    def dimensions(self):
        """Return image dimensions as a tuple."""
        dim = self.get_dimensions()
        return tuple(dim[1:1 + dim[0]]) # dim[0] tells us how many dimensions

    def get_dimensions(self):
        """Get dimensions.

        Returns nifti format.
        i.e., dim[0]=ndims, dim[1]=nx, dim[2]=ny,...
        """
        dim = numpy.ndarray((8,), dtype=cpp_int_dtype())
        self.handle.cReg_NiftiImageData_get_dimensions(dim)
        return dim

    def get_voxel_sizes(self):
        """Get voxel sizes.

        Returns nifti format.
        i.e., dim[0]=?, dim[1]=dx, dim[2]=dy,...
        """
        out = numpy.ndarray((8,), dtype=numpy.float32)
        self.handle.cReg_NiftiImageData_get_voxel_sizes(out)
        return out

    def dot(self, other):
        '''
        Returns the dot product of the container data with another container
        data or numpy array viewed as vectors.
        other: NiftiImageData or numpy array.
        '''
        if not (issubclass(type(other), type(self))):
            self_copy = self.clone()
            self_copy.fill(other)
            other = self_copy
        return super().dot(other)

    def add(self, other, out=None):
        '''
        Addition for NiftiImageData containers.

        If other is a NiftiData or numpy array, returns the sum of data
        stored in self and other viewed as vectors.
        If other is a scalar, returns the same with the second vector filled
        with the value of other.
        other: NiftiImageData or numpy array or scalar.
        out:   NiftiImageData to store the result to.
        '''
        if not (issubclass(type(other), type(self)) or isinstance(other, (Number, numpy.number))):
            self_copy = self.clone()
            self_copy.fill(other)
            other = self_copy
        return super().add(other, out)

    def subtract(self, other, out=None):
        '''
        Subtraction for NiftiImageData containers.

        If other is a NiftiData or numpy array, returns the sum of data
        stored in self and other viewed as vectors.
        If other is a scalar, returns the same with the second vector filled
        with the value of other.
        other: NiftiImageData or numpy array or scalar.
        out:   NiftiImageData to store the result to.
        '''
        if not (issubclass(type(other), type(self)) or isinstance(other, (Number, numpy.number))):
            self_copy = self.clone()
            self_copy.fill(other)
            other = self_copy
        return super().subtract(other, out)

    def binary(self, other, f, out=None):
        '''Applies function f(x,y) element-wise to self and other.

        other: NiftiImageData or numpy array or Number
        f: the name of the function to apply, Python str.
        '''
        if not (issubclass(type(other), type(self)) or isinstance(other, (Number, numpy.number))):
            self_copy = self.clone()
            self_copy.fill(other)
            other = self_copy
        return super().binary(other, f, out)

    def fill(self, val):
        """Fill image with single value or numpy array."""
        if isinstance(val, numpy.ndarray):
            if val.dtype is numpy.dtype('float32'):
                v = val
            else:
                v = val.astype(numpy.float32)
            if not v.flags['F_CONTIGUOUS']:
                v = numpy.asfortranarray(v)
            self.handle.cReg_NiftiImageData_fill_arr(v)
        elif isinstance(val, (float, int)):
            self.handle.cReg_NiftiImageData_fill(float(val))
        else:
            raise error("""wrong fill value. Should be numpy.ndarray,
                        float or int""")

    def deep_copy(self):
        """Deep copy image."""
        if self.name == 'NiftiImageData':
            image = NiftiImageData()
        elif self.name == 'NiftiImageData3D':
            image = NiftiImageData3D()
        elif self.name == 'NiftiImageData3DTensor':
            image = NiftiImageData3DTensor()
        elif self.name == 'NiftiImageData3DDeformation':
            image = NiftiImageData3DDeformation()
        elif self.name == 'NiftiImageData3DDisplacement':
            image = NiftiImageData3DDisplacement()
        image.handle.cReg_NiftiImageData_deep_copy(self)
        return image

    def as_array(self):
        """Get data as numpy array."""
        dim = self.get_dimensions()
        dim = dim[1:dim[0] + 1]
        array = numpy.ndarray(dim, dtype=numpy.float32, order='F')
        self.handle.cReg_NiftiImageData_as_array(array)
        return numpy.ascontiguousarray(array)

    def get_original_datatype(self):
        """Get original image datatype.

        (internally everything is converted to float).
        """
        datatype = self.handle.cReg_NiftiImageData_get_original_datatype()
        return int(datatype)

    def crop(self, min_, max_):
        """Crop image.

        Give minimum and maximum indices.
        Min and max indicies can be anywhere between
        (x,y,z) and (x,y,z,t,u,v,w).
        Use values of -1 for no change.
        """
        if len(min_) < 3 or len(min_) > 7:
            raise AssertionError("""Min bounds should be at least (x,y,z),
                                 and up to (x,y,z,t,u,v,w)""")
        if len(max_) < 3 or len(max_) > 7:
            raise AssertionError("""Max bounds should be at least (x,y,z),
                                 and up to (x,y,z,t,u,v,w)""")
        # Fill in any missing indices with -1's
        min_.extend([-1] * (7 - len(min_)))
        max_.extend([-1] * (7 - len(max_)))
        min_np = numpy.array(min_, dtype=cpp_int_dtype())
        max_np = numpy.array(max_, dtype=cpp_int_dtype())
        self.handle.cReg_NiftiImageData_crop(min_np, max_np)

    def pad(self, min_, max_, val=0):
        """Pad image.

        Give minimum and maximum indices.
        Min and max indicies can be anywhere between
        (x,y,z) and (x,y,z,t,u,v,w).
        Use values of -1 for no change.
        """
        if len(min_) < 3 or len(min_) > 7:
            raise AssertionError("""Min bounds should be at least (x,y,z),
                                 and up to (x,y,z,t,u,v,w)""")
        if len(max_) < 3 or len(max_) > 7:
            raise AssertionError("""Max bounds should be at least (x,y,z),
                                 and up to (x,y,z,t,u,v,w)""")
        # Fill in any missing indices with -1's
        min_.extend([-1] * (7 - len(min_)))
        max_.extend([-1] * (7 - len(max_)))
        min_np = numpy.array(min_, dtype=cpp_int_dtype())
        max_np = numpy.array(max_, dtype=cpp_int_dtype())
        self.handle.cReg_NiftiImageData_pad(min_np, max_np, float(val))

    def print_header(self):
        """Print nifti header metadata."""
        print(self.get_header())

    def get_header(self):
        """Print nifti header metadata."""
        vec = SIRF.DataHandleVector()
        vec.push_back(self.handle)
        #vec.handle.cReg_NiftiImageData_print_headers()
        s = vec.handle.cReg_NiftiImageData_get_headers()
        return str(s)

    def same_object(self):
        """See DataContainer method."""
        obj = ImageData()
        obj.handle = None
        return obj
        #return NiftiImageData()

    def set_voxel_spacing(self, spacing, interpolation_order):
        """Set the voxel spacing.

        Requires resampling image,
        and so interpolation order is required.
        As per NiftyReg, interpolation_order can be either 0, 1 or 3
        meaning nearest neighbor, linear or cubic spline interpolation.
        """
        if len(spacing) != 3:
            raise AssertionError("New spacing should be array of 3 numbers.")
        self.handle.cReg_NiftiImageData_set_voxel_spacing(float(spacing[0]), float(spacing[1]), float(spacing[2]),
                                                          int(interpolation_order))

    def get_contains_nans(self):
        """Return true if image contains any voxels with NaNs."""
        return parms.bool_par(self.handle, 'NiftiImageData', 'contains_nans')

    def normalise_zero_and_one(self):
        """Normalise image between 0 and 1."""
        self.handle.cReg_NiftiImageData_normalise_zero_and_one()

    def standardise(self):
        """Standardise (subtract mean and divide by standard deviation)."""
        self.handle.cReg_NiftiImageData_standardise()

    def get_inner_product(self, other):
        """Get inner product between two images. Must be same size."""
        if not isinstance(other, NiftiImageData):
            raise AssertionError()
        inner_product = self.handle.cReg_NiftiImageData_get_inner_product(other)
        return float(inner_product)

    def kernel_convolution(self, sigma, convolution_type=2):
        """Kernel convolution.

        convolution_type:
            - MEAN_KERNEL         = 0
            - LINEAR_KERNEL       = 1
            - GAUSSIAN_KERNEL     = 2
            - CUBIC_SPLINE_KERNEL = 3
        """
        self.handle.cReg_NiftiImageData_kernel_convolution(float(sigma), int(convolution_type))

    @staticmethod
    def print_headers(to_print):
        """Print nifti header metadata of one or multiple nifti images."""
        print(NiftiImageData.get_headers(to_print))

    @staticmethod
    def get_headers(to_print):
        if not all(isinstance(n, NiftiImageData) for n in to_print):
            raise AssertionError()
        vec = SIRF.DataHandleVector()
        for n in to_print:
            vec.push_back(n.handle)
        #vec.handle.cReg_NiftiImageData_print_headers()
        s = vec.handle.cReg_NiftiImageData_get_headers()
        return str(s)

    @staticmethod
    def construct_from_complex_image(complex_im):
        """Construct two NiftiImageData from a complex image."""
        if not isinstance(complex_im, SIRF.ImageData):
            raise AssertionError()
        im_real = NiftiImageData()
        im_imag = NiftiImageData()
        im_real.handle = complex_im.handle.cReg_NiftiImageData_from_complex_ImageData_real_component()
        im_imag.handle = complex_im.handle.cReg_NiftiImageData_from_complex_ImageData_imag_component()
        return [im_real, im_imag]

    @staticmethod
    def are_equal_to_given_accuracy(im1, im2, accuracy):
        """Check if two images match to a given accuracy."""
        if not isinstance(im1, NiftiImageData) or not \
                isinstance(im2, NiftiImageData):
            raise AssertionError()
        if not im1.handle.valid or not im2.handle.valid:
            error('Cannot compare images as at least one is uninitialised')
        h = im1.handle.cReg_NiftiImageData_are_equal_to_given_accuracy(im2, float(accuracy))
        return int(h)

    @property
    def shape(self):
        return self.dimensions()

    @property
    def __array_interface__(self):
        """As per https://numpy.org/doc/stable/reference/arrays.interface.html"""
        if not self.supports_array_view:
            raise ContiguousError("please make an array-copy first with `asarray(copy=True)` or `as_array()`")
        # NiftiImageData stores voxels values in a 3D Fortran-style array,
        # hence strides need to be arranged accordingly (default is C-style).
        shape = self.shape
        strides = ()
        stride = 4
        for d in shape:
            strides += (stride,)
            stride *= d
        return {
            'shape': shape, 'typestr': '<f4', 'strides': strides, 'version': 3,
            'data': (parms.size_t_par(self.handle, 'NiftiImageData', 'address'), False)}


class NiftiImageData3D(NiftiImageData):
    """Class for 3D nifti image data."""
    def __init__(self, src=None):
        """init."""
        self.handle = None
        self.name = 'NiftiImageData3D'
        if src is None:
            self.handle = Handle(None, -1).cReg_newObject(self.name)
        elif isinstance(src, str):
            self.handle = Handle(None, -1).cReg_objectFromFile(self.name, src)
        elif isinstance(src, SIRF.ImageData):
            # src is ImageData
            self.handle = src.handle.cReg_NiftiImageData_from_SIRFImageData(1)
        else:
            raise error('Wrong source in NiftiImageData3D constructor')


class ImageData(NiftiImageData3D):
    """Alias class for nifti image data."""


class NiftiImageData3DTensor(NiftiImageData):
    """Class for 3D tensor nifti image data."""
    def __init__(self, src1=None, src2=None, src3=None):
        """init."""
        self.handle = None
        self.name = 'NiftiImageData3DTensor'
        if src1 is None:
            self.handle = Handle(None, -1).cReg_newObject(self.name)
        elif isinstance(src1, str):
            self.handle = Handle(None, -1).cReg_objectFromFile(self.name, src1)
        elif all(isinstance(i, NiftiImageData3D) for i in (src1, src2, src3)):
            self.handle = Handle(None, -1).cReg_NiftiImageData3DTensor_construct_from_3_components(
                self.name, src1, src2, src3)
        else:
            raise error('Wrong source in NiftiImageData3DTensor constructor')

    def write_split_xyz_components(self, filename, datatype=-1):
        """Save to file.

        See nifti1.h for datatypes
        (e.g., float (NIFTI_TYPE_FLOAT32) = 16).
        Image's original datatpye is used by default.
        """
        if not isinstance(filename, str):
            raise AssertionError()
        self.handle.cReg_NiftiImageData3DTensor_write_split_xyz_components(filename, datatype)

    def create_from_3D_image(self, src):
        """Create tensor/deformation/displacement field from 3D image."""
        if not isinstance(src, NiftiImageData3D):
            raise AssertionError()
        self.handle.cReg_NiftiImageData3DTensor_create_from_3D_image(src)

    def flip_component(self, dim):
        """Flip component of nu."""
        if 0 < dim or dim > 2:
            raise AssertionError("Dimension to flip should be between 0 and 2.")
        self.handle.cReg_NiftiImageData3DTensor_flip_component(dim)


class NiftiImageData3DDisplacement(NiftiImageData3DTensor, _Transformation):
    """Class for 3D displacement nifti image data.

    Displacement fields (as opposed to Deformation fields) describe the change
    (in real world units) of the pixel locations between images.
    A displacement field of an identity transformation will be of zero value.
    """
    def __init__(self, src1=None, src2=None, src3=None):
        """init."""
        self.handle = None
        self.name = 'NiftiImageData3DDisplacement'
        if src1 is None:
            self.handle = Handle(None, -1).cReg_newObject(self.name)
        elif isinstance(src1, str):
            self.handle = Handle(None, -1).cReg_objectFromFile(self.name, src1)
        elif all(isinstance(i, NiftiImageData3D) for i in (src1, src2, src3)):
            self.handle = Handle(None, -1).cReg_NiftiImageData3DTensor_construct_from_3_components(
                self.name, src1, src2, src3)
        elif isinstance(src1, NiftiImageData3DDeformation):
            self.handle = src1.handle.cReg_NiftiImageData3DDisplacement_create_from_def()
        else:
            raise error('Wrong source in NiftiImageData3DDisplacement constructor')


class NiftiImageData3DDeformation(NiftiImageData3DTensor, _Transformation):
    """Class for 3D deformation nifti image data.

    Deformation fields (as opposed to Displacement fields) describe the
    absolute position (in real world units) of the pixel locations on the
    reference image.
    A deformation field of an identity transformation will contain the
    location of each of the pixels centroids in the world coordinates.
    """
    def __init__(self, src1=None, src2=None, src3=None):
        """init."""
        self.handle = None
        self.name = 'NiftiImageData3DDeformation'
        if src1 is None:
            self.handle = Handle(None, -1).cReg_newObject(self.name)
        elif isinstance(src1, str):
            self.handle = Handle(None, -1).cReg_objectFromFile(self.name, src1)
        elif all(isinstance(i, NiftiImageData3D) for i in (src1, src2, src3)):
            self.handle = Handle(None, -1).cReg_NiftiImageData3DTensor_construct_from_3_components(
                self.name, src1, src2, src3)
        elif isinstance(src1, NiftiImageData3DDisplacement):
            self.handle = src1.handle.cReg_NiftiImageData3DDeformation_create_from_disp()
        else:
            raise error('Wrong source in NiftiImageData3DDeformation constructor')

    def get_inverse(self, floating=None):
        """
        Get inverse (potentially based on another image).

        Why would you want to base it on another image? Well, we might have a
        deformation that takes us from image A to B. We'll probably want the
        inverse to take us from image B back to A. In this case, use
        get_inverse(A). This is because the the deformation field is defined
        for the reference image. In the second case, A is the reference,
        and B is the floating image.
        """
        if floating is None:
            floating = self
        if not isinstance(floating, SIRF.ImageData):
            raise AssertionError()
        if self is None:
            raise AssertionError()
        output = NiftiImageData3DDeformation()
        output.handle = self.handle.cReg_NiftiImageData3DDeformation_get_inverse(floating)
        return output

    @staticmethod
    def compose_single_deformation(trans, ref):
        """Compose transformations into single deformation."""
        if not isinstance(ref, NiftiImageData3D):
            raise AssertionError()
        if not all(isinstance(n, _Transformation) for n in trans):
            raise AssertionError()
        if len(trans) == 1:
            return trans[0].get_as_deformation_field(ref)
        # This is ugly. Store each type in a single string (need to do this
        # because I can't get virtual methods to work for multiple inheritance
        # (deformation/displacement are both nifti images and transformations).
        types = ''
        for n in trans:
            if isinstance(n, AffineTransformation):
                types += '1'
            elif isinstance(n, NiftiImageData3DDisplacement):
                types += '2'
            elif isinstance(n, NiftiImageData3DDeformation):
                types += '3'
        # Convert transformations into SIRF vector
        vec = SIRF.DataHandleVector()
        for n in trans:
            vec.push_back(n.handle)
        z = NiftiImageData3DDeformation()
        z.handle = ref.handle.cReg_NiftiImageData3DDeformation_compose_single_deformation(types, vec)
        return z


class _Registration(ABC):
    """Abstract base class for registration."""
    def __init__(self):
        """init."""
        self.handle = None
        self.name = 'Registration'
        self.reference_image = None

    def set_reference_image(self, reference_image):
        """Set the reference image."""
        if not isinstance(reference_image, SIRF.ImageData):
            raise AssertionError()
        self.reference_image = reference_image
        parms.set_parameter(self.handle, 'Registration', 'reference_image', reference_image.handle)

    def set_floating_image(self, floating_image):
        """Set the floating image. Will clear any previous floating images."""
        if not isinstance(floating_image, SIRF.ImageData):
            raise AssertionError()
        parms.set_parameter(self.handle, 'Registration', 'floating_image', floating_image.handle)

    def add_floating_image(self, floating_image):
        """Add floating image."""
        if not isinstance(floating_image, SIRF.ImageData):
            raise AssertionError()
        self.handle.cReg_Registration_add_floating(floating_image)

    def set_reference_image_filename(self, filename):
        """Set reference image filename."""
        if not isinstance(filename, str):
            raise AssertionError()
        self.reference_image = NiftiImageData(filename)
        self.handle.cReg_Registration_set_reference_image_filename(filename)

    def set_floating_image_filename(self, filename):
        """Set floating image filename."""
        if not isinstance(filename, str):
            raise AssertionError()
        self.handle.cReg_Registration_set_floating_image_filename(filename)

    def add_floating_image_filename(self, filename):
        """Add floating image filename."""
        if not isinstance(filename, str):
            raise AssertionError()
        self.handle.cReg_Registration_add_floating_image_filename(filename)

    def clear_floating_images(self):
        """Clear floating images."""
        self.handle.cReg_Registration_clear_floatings()

    def get_output(self, idx=0):
        """Get the registered image."""
        output = self.reference_image.same_object()
        output.handle = self.handle.cReg_Registration_get_output(int(idx))
        return output

    def process(self):
        """Run the registration."""
        self.handle.cReg_Registration_process()

    def get_deformation_field_forward(self, idx=0):
        """Get the forward deformation field image."""
        output = NiftiImageData3DDeformation()
        output.handle = self.handle.cReg_Registration_get_deformation_displacement_image(
            'forward_deformation', int(idx))
        return output

    def get_deformation_field_inverse(self, idx=0):
        """Get the inverse deformation field image."""
        output = NiftiImageData3DDeformation()
        output.handle = self.handle.cReg_Registration_get_deformation_displacement_image(
            'inverse_deformation', int(idx))
        return output

    def get_displacement_field_forward(self, idx=0):
        """Get the forward displacement field image."""
        output = NiftiImageData3DDisplacement()
        output.handle = self.handle.cReg_Registration_get_deformation_displacement_image(
            'forward_displacement', int(idx))
        return output

    def get_displacement_field_inverse(self, idx=0):
        """Get the inverse displacement field image."""
        output = NiftiImageData3DDisplacement()
        output.handle = self.handle.cReg_Registration_get_deformation_displacement_image(
            'inverse_displacement', int(idx))
        return output


class _NiftyRegistration(_Registration):
    """Abstract base class for NiftyReg registration."""
    def __init__(self):
        """init."""
        super().__init__()
        self.name = 'NiftyRegistration'

    def set_parameter_file(self, filename):
        """Set the parameter filename."""
        parms.set_char_par(self.handle, 'NiftyRegistration', 'parameter_file', filename)

    def set_reference_mask(self, reference_mask):
        """Set the reference mask."""
        if not isinstance(reference_mask, SIRF.ImageData):
            raise AssertionError()
        parms.set_parameter(self.handle, 'NiftyRegistration', 'reference_mask', reference_mask.handle)

    def set_floating_mask(self, floating_mask):
        """Set the floating mask."""
        if not isinstance(floating_mask, SIRF.ImageData):
            raise AssertionError()
        parms.set_parameter(self.handle, 'NiftyRegistration', 'floating_mask', floating_mask.handle)

    def set_parameter(self, par, arg1="", arg2=""):
        """Set string parameter.

        Check if any set methods match the method
        given by par. If so, set the value given by arg. Convert to float/int
        etc., as necessary. Up to 2 arguments, leave blank if unneeded.
        These are applied after parsing the parameter file.
        """
        self.handle.cReg_NiftyRegistration_set_parameter(par, arg1, arg2)


class NiftyAladinSym(_NiftyRegistration):
    """Registration using NiftyReg aladin."""
    def __init__(self):
        """init."""
        super().__init__()
        self.name = 'NiftyAladinSym'
        self.handle = Handle(None, -1).cReg_newObject(self.name)

    def get_transformation_matrix_forward(self):
        """Get forward transformation matrix."""
        tm = AffineTransformation()
        tm.handle = self.handle.cReg_NiftyAladin_get_TM('forward')
        return tm

    def get_transformation_matrix_inverse(self):
        """Get inverse transformation matrix."""
        tm = AffineTransformation()
        tm.handle = self.handle.cReg_NiftyAladin_get_TM('inverse')
        return tm

    @staticmethod
    def print_all_wrapped_methods():
        """Print all wrapped methods."""
        print("""In C++, this class is templated. \"dataType\"
              corresponds to \"float\" for Matlab and python.""")
        Handle(None, -1).cReg_NiftyRegistration_print_all_wrapped_methods('NiftyAladinSym')


class NiftyF3dSym(_NiftyRegistration):
    """Registration using NiftyReg f3d."""
    def __init__(self):
        """init."""
        super().__init__()
        self.name = 'NiftyF3dSym'
        self.handle = Handle(None, -1).cReg_newObject(self.name)

    def set_floating_time_point(self, floating_time_point):
        """Set floating time point."""
        parms.set_int_par(self.handle, self.name, 'floating_time_point', floating_time_point)

    def set_reference_time_point(self, reference_time_point):
        """Set reference time point."""
        parms.set_int_par(self.handle, self.name, 'reference_time_point', reference_time_point)

    def set_initial_affine_transformation(self, src):
        """Set initial affine transformation."""
        if not isinstance(src, AffineTransformation):
            raise AssertionError()
        parms.set_parameter(self.handle, self.name, 'initial_affine_transformation', src.handle)

    @staticmethod
    def print_all_wrapped_methods():
        """Print all wrapped methods."""
        print("""In C++, this class is templated. \"dataType\"
              corresponds to \"float\" for Matlab and python.""")
        Handle(None, -1).cReg_NiftyRegistration_print_all_wrapped_methods('NiftyF3dSym')


if SIRF_HAS_SPM:

    class SPMRegistration(_Registration):
        """Registration using SPM."""
        def __init__(self):
            """init."""
            super().__init__()
            self.name = 'SPMRegistration'
            self.handle = Handle(None, -1).cReg_newObject(self.name)

        def get_transformation_matrix_forward(self, idx=0):
            """Get forward transformation matrix."""
            tm = AffineTransformation()
            tm.handle = self.handle.cReg_SPMRegistration_get_TM('forward', int(idx))
            return tm

        def get_transformation_matrix_inverse(self, idx=0):
            """Get inverse transformation matrix."""
            tm = AffineTransformation()
            tm.handle = self.handle.cReg_SPMRegistration_get_TM('inverse', int(idx))
            return tm

        def set_working_folder(self, working_folder):
            """Set working folder."""
            parms.set_char_par(self.handle, self.name, 'working_folder', working_folder)

        def set_working_folder_file_overwrite(self, working_folder_file_overwrite=True):
            """Set file overwrite in working folder."""
            if working_folder_file_overwrite:
                working_folder_file_overwrite = 1
            else:
                working_folder_file_overwrite = 0
            parms.set_int_par(self.handle, self.name, 'working_folder_file_overwrite', working_folder_file_overwrite)

        def set_delete_temp_files(self, delete_temp_files=True):
            """Delete temporary files."""
            if delete_temp_files:
                delete_temp_files = 1
            else:
                delete_temp_files = 0
            parms.set_int_par(self.handle, self.name, 'delete_temp_files', delete_temp_files)


class NiftyResampler:
    """Resample using NiftyReg."""
    def __init__(self):
        """init."""
        self.name = 'NiftyResampler'
        self.handle = Handle(None, -1).cReg_newObject(self.name)
        self.reference_image = None
        self.floating_image = None

    def set_reference_image(self, reference_image):
        """Set reference image.

        This is the image that would be the reference
        if you were doing a forward transformation.
        """
        if not isinstance(reference_image, SIRF.ImageData):
            raise AssertionError()
        self.reference_image = reference_image
        parms.set_parameter(self.handle, self.name, 'reference_image', reference_image.handle)

    def set_floating_image(self, floating_image):
        """Set floating image.

        This is the image that would be the floating if
        you were doing a forward transformation.
        """
        if not isinstance(floating_image, SIRF.ImageData):
            raise AssertionError()
        self.floating_image = floating_image
        parms.set_parameter(self.handle, self.name, 'floating_image', floating_image.handle)

    def add_transformation(self, src):
        """Add transformation."""
        if isinstance(src, AffineTransformation):
            self.handle.cReg_NiftyResampler_add_transformation(src, 'affine')
        elif isinstance(src, NiftiImageData3DDisplacement):
            self.handle.cReg_NiftyResampler_add_transformation(src, 'displacement')
        elif isinstance(src, NiftiImageData3DDeformation):
            self.handle.cReg_NiftyResampler_add_transformation(src, 'deformation')
        else:
            raise AssertionError()

    def clear_transformations(self):
        """Clear transformations."""
        if self.handle.valid:
            self.handle.cReg_NiftyResampler_clear_transformations()

    def set_interpolation_type(self, interp_type):
        """Set interpolation type.

        0=nearest neighbour, 1=linear, 3=cubic, 4=sinc.
        """
        if not isinstance(interp_type, int):
            raise AssertionError()
        parms.set_int_par(self.handle, self.name, 'interpolation_type', interp_type)

    def set_interpolation_type_to_nearest_neighbour(self):
        """Set interpolation type to nearest neighbour."""
        parms.set_int_par(self.handle, self.name, 'interpolation_type', 0)

    def set_interpolation_type_to_linear(self):
        """Set interpolation type to linear."""
        parms.set_int_par(self.handle, self.name, 'interpolation_type', 1)

    def set_interpolation_type_to_cubic_spline(self):
        """Set interpolation type to cubic spline."""
        parms.set_int_par(self.handle, self.name, 'interpolation_type', 3)

    def set_interpolation_type_to_sinc(self):
        """Set interpolation type to sinc."""
        parms.set_int_par(self.handle, self.name, 'interpolation_type', 4)

    def set_padding_value(self, val):
        """Set padding value."""
        parms.set_float_par(self.handle, self.name, 'padding', val)

    def norm(self, num_iter=2, verb=0):
        '''Computes the norm of the forward projection operator.
        '''
        # reference and floating images need to be real.
        if self.reference_image.is_complex():
            reference_image, _ = NiftiImageData.construct_from_complex_image(self.reference_image)
            parms.set_parameter(self.handle, self.name, 'reference_image', reference_image.handle)
        if self.floating_image.is_complex():
            floating_image, _ = NiftiImageData.construct_from_complex_image(self.floating_image)
            parms.set_parameter(self.handle, self.name, 'floating_image', floating_image.handle)
        handle = self.handle.cReg_NiftyResampler_norm(num_iter, verb)
        r = float(handle)
        # Restore reference and floating images.
        if self.reference_image.is_complex():
            parms.set_parameter(self.handle, self.name, 'reference_image', self.reference_image.handle)
        if self.floating_image.is_complex():
            parms.set_parameter(self.handle, self.name, 'floating_image', self.floating_image.handle)
        return r

    def process(self):
        """Process.

        Equivalent of calling forward(floating_image).
        Use get_output to get resampled image.
        """
        self.handle.cReg_NiftyResampler_process()

    def get_output(self):
        """Get output."""
        image = self.reference_image.same_object()
        image.handle = parms.parameter_handle(self.handle, self.name, 'output')
        return image

    def forward(self, x, out=None):
        """
        Forward transformation.

        Usage:
            output = forward(x), OR
            forward(x=input,out=output)
        """
        if out is None:
            if self.reference_image is None:
                raise RuntimeError('reference_image needs to be set first, or you should specify an out argument')
            out = self.reference_image.clone()
        # Check image validity
        if not isinstance(x, SIRF.ImageData):
            raise TypeError('{} expecting input as SIRF.ImageData, got {}'.format(self.__class__.__name__, type(x)))

        if not isinstance(out, SIRF.ImageData):
            raise TypeError('{} expecting output as SIRF.ImageData, got {}'.format(self.__class__.__name__, type(out)))
        # Forward
        out.handle.cReg_NiftyResampler_forward(x, self)
        return out

    def direct(self, x, out=None):
        """Alias to forward."""
        return self.forward(x=x, out=out)

    def adjoint(self, x, out=None):
        """
        Adjoint transformation.

        Usage:
            output = adjoint(x), OR
            adjoint(x=input,out=output)
        """
        if out is None:
            out = self.floating_image.clone()
        if not isinstance(x, SIRF.ImageData):
            raise TypeError("""{} expecting input as subclass of
                            SIRF.ImageData, got {}""".format(self.__class__.__name__, type(x)))

        if not isinstance(out, SIRF.ImageData):
            raise TypeError('{} expecting output as SIRF.ImageData, got {}'.format(self.__class__.__name__, type(out)))
        # Forward
        out.handle.cReg_NiftyResampler_adjoint(x, self)
        return out

    def backward(self, x, out=None):
        """
        Backward transformation.

        Alias of adjoint to align terms with
        AcquisitionModel's forward and backward.
        """
        return self.adjoint(x=x, out=out)

    def is_linear(self):
        """Return whether the transformation is linear."""
        return True

    def domain_geometry(self):
        """Get domain geometry."""
        return self.floating_image

    def range_geometry(self):
        """Get range geometry."""
        return self.reference_image


"""This alias for the NiftyResampler is deprecated."""
NiftyResample = NiftyResampler


class ImageWeightedMean:
    """Class for performing weighted mean of images."""
    def __init__(self):
        """init."""
        self.name = 'ImageWeightedMean'
        self.handle = Handle(None, -1).cReg_newObject(self.name)

    def add_image(self, image, weight):
        """Add an image and its corresponding weight.

        Image should be via filename or NiftiImageData.
        """
        if isinstance(image, NiftiImageData):
            self.handle.cReg_ImageWeightedMean_add_image(image, weight)
        elif isinstance(image, str):
            self.handle.cReg_ImageWeightedMean_add_image_filename(image, weight)
        else:
            raise error("""sirf.Reg.ImageWeightedMean.add_image: image must be
                        NiftiImageData or filename.""")

    def process(self):
        """Process."""
        self.handle.cReg_ImageWeightedMean_process()

    def get_output(self):
        """Get output."""
        image = NiftiImageData()
        image.handle = parms.parameter_handle(self.handle, self.name, 'output')
        return image


class AffineTransformation(_Transformation):
    """Class for affine transformations."""
    def __init__(self, src1=None, src2=None):
        """init."""
        self.handle = None
        self.name = 'AffineTransformation'
        if src1 is None:
            self.handle = Handle(None, -1).cReg_newObject(self.name)
        elif isinstance(src1, str):
            self.handle = Handle(None, -1).cReg_objectFromFile(self.name, src1)
        elif isinstance(src1, numpy.ndarray) and src2 is None:
            src1 = format_numpy_array_for_setter(src1)
            if src1.shape != (4, 4):
                raise AssertionError()
            # Need to transpose relative to MATLAB
            trans = numpy.zeros((4, 4), dtype=numpy.float32)
            for i in range(4):
                for j in range(4):
                    trans[i, j] = src1[j, i]
                self.handle = Handle(None, -1).cReg_AffineTransformation_construct_from_TM(trans)
        elif isinstance(src1, numpy.ndarray) and src2 is not None and \
                isinstance(src2, Quaternion):
            src1 = format_numpy_array_for_setter(src1)
            self.handle = src2.handle.cReg_AffineTransformation_construct_from_trans_and_quaternion(src1)
        elif isinstance(src1, numpy.ndarray) and isinstance(src2, numpy.ndarray):
            src1 = format_numpy_array_for_setter(src1)
            src2 = format_numpy_array_for_setter(src2)
            self.handle = Handle(None, -1).cReg_AffineTransformation_construct_from_trans_and_euler(src1, src2)
        else:
            raise error("""AffineTransformation accepts no args, filename,
                        4x4 array or translation with quaternion.""")

    def __eq__(self, other):
        """Overload comparison operator."""
        if not isinstance(other, AffineTransformation):
            raise AssertionError()
        h = self.handle.cReg_AffineTransformation_equal(other)
        return int(h)

    def __ne__(self, other):
        """Overload comparison operator."""
        return not self == other

    def __mul__(self, other):
        """Overload multiplication operator."""
        if not isinstance(other, AffineTransformation):
            raise AssertionError()
        mat = AffineTransformation()
        mat.handle = self.handle.cReg_AffineTransformation_mul(other)
        return mat

    def deep_copy(self):
        """Deep copy."""
        mat = AffineTransformation()
        mat.handle = self.handle.cReg_AffineTransformation_deep_copy()
        return mat

    def write(self, filename):
        """Save to file."""
        self.handle.cReg_AffineTransformation_write(filename)

    def get_determinant(self):
        """Get determinant."""
        return parms.float_par(self.handle, self.name, 'determinant')

    def as_array(self):
        """Get forward transformation matrix."""
        tm = numpy.ndarray((4, 4), dtype=numpy.float32)
        self.handle.cReg_AffineTransformation_as_array(tm)
        return tm

    def get_inverse(self):
        """Get inverse matrix."""
        tm = AffineTransformation()
        tm.handle = self.handle.cReg_AffineTransformation_get_inverse()
        return tm

    def get_Euler_angles(self):
        """Get Euler angles of transformation matrix (XYZ)."""
        eul = numpy.ndarray(3, dtype=numpy.float32)
        self.handle.cReg_AffineTransformation_get_Euler_angles(eul)
        return eul

    def get_quaternion(self):
        """Get quaternion."""
        quat_zeros = numpy.array([0., 0., 0., 0.], dtype=numpy.float32)
        quat = Quaternion(quat_zeros)
        quat.handle = self.handle.cReg_AffineTransformation_get_quaternion()
        return quat

    @staticmethod
    def get_identity():
        """Get identity matrix."""
        mat = AffineTransformation()
        mat.handle = Handle(None, -1).cReg_AffineTransformation_get_identity()
        return mat

    @staticmethod
    def get_average(to_average):
        """Get average of transformations."""
        if not all(isinstance(n, AffineTransformation) for n in to_average):
            raise AssertionError("""AffineTransformation:get_average() input
                                 list should only contain
                                 AffineTransformations.""")
        tm = AffineTransformation()
        vec = SIRF.DataHandleVector()
        for n in to_average:
            vec.push_back(n.handle)
        tm.handle = vec.handle.cReg_AffineTransformation_get_average()
        return tm


class Quaternion:
    """Class for quaternions."""
    def __init__(self, src=None):
        """init."""
        self.handle = None
        self.name = 'Quaternion'
        if isinstance(src, numpy.ndarray):
            if src.size != 4:
                raise AssertionError("""Quaternion constructor from numpy
                                     array is wrong size.""")
            if src.dtype is not numpy.float32:
                src = src.astype(numpy.float32)
            self.handle = Handle(None, -1).cReg_Quaternion_construct_from_array(src)
        elif isinstance(src, AffineTransformation):
            self.handle = src.handle.cReg_Quaternion_construct_from_AffineTransformation()
        else:
            raise error('Wrong source in quaternion constructor')

    def as_array(self):
        """Get quaternion as array."""
        arr = numpy.ndarray(4, dtype=numpy.float32)
        self.handle.cReg_Quaternion_as_array(arr)
        return arr

    @staticmethod
    def get_average(to_average):
        """Get average of quaternions."""
        if not all(isinstance(n, Quaternion) for n in to_average):
            raise AssertionError()
        quat_zeros = numpy.array([0., 0., 0., 0.], dtype=numpy.float32)
        quat = Quaternion(quat_zeros)
        if not all(isinstance(n, Quaternion) for n in to_average):
            raise AssertionError()
        vec = SIRF.DataHandleVector()
        for n in to_average:
            vec.push_back(n.handle)
        quat.handle = vec.handle.cReg_Quaternion_get_average()
        return quat
