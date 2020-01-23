"""
Object-Oriented wrap for the cReg-to-Python interface pyreg.py
"""

# CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
# Copyright 2018 - 2019 University College London
#
# This is software developed for the Collaborative Computational
# Project in Positron Emission Tomography and Magnetic Resonance imaging
# (http://www.ccppetmr.ac.uk/).
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
import time
import numbers
import deprecation

from pUtilities import *
from sirf import SIRF
import pyiutilities as pyiutil
import pyreg

import sirf.Reg_params as parms


try:
    input_ = raw_input
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

class MessageRedirector(object):
    """
    Class for registration printing redirection to files/stdout/stderr.
    """
    def __init__(self, info=None, warn='stdout', errr='stdout'):
        """
        Creates MessageRedirector object that redirects Reg's ouput
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
        if self.info_case == 0:
            try_calling(pyreg.deleteTextPrinter(self.info))
        else:
            try_calling(pyreg.deleteTextWriter(self.info))
        pyreg.closeChannel(0, self.info)
        if self.warn_case == 0:
            try_calling(pyreg.deleteTextPrinter(self.warn))
        else:
            try_calling(pyreg.deleteTextWriter(self.warn))
        pyreg.closeChannel(1, self.warn)
        if self.errr_case == 0:
            try_calling(pyreg.deleteTextPrinter(self.errr))
        else:
            try_calling(pyreg.deleteTextWriter(self.errr))
        pyreg.closeChannel(2, self.errr)
###########################################################


class _Transformation(ABC):
    """
    Abstract base class for transformations.
    """
    def __init__(self):
        self.handle = None
        self.name = 'Transformation'

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def get_as_deformation_field(self, ref):
        """Get any type of transformation as a deformation field.
        This is useful for joining them together. Require a reference
        image for converting transformation matrices to deformations."""
        if self.handle is None:
            raise AssertionError()
        if not isinstance(ref, NiftiImageData3D):
            raise AssertionError()
        output = NiftiImageData3DDeformation()
        output.handle = pyreg.cReg_Transformation_get_as_deformation_field(self.handle, self.name, ref.handle)
        check_status(output.handle)
        return output


class NiftiImageData(SIRF.ImageData):
    """
    General class for nifti image data.
    """
    def __init__(self, src=None):
        self.handle = None
        self.name = 'NiftiImageData'
        if src is None:
            self.handle = pyreg.cReg_newObject(self.name)
        elif isinstance(src, str):
            self.handle = pyreg.cReg_objectFromFile(self.name, src)
        else:
            raise error('Wrong source in NiftiImageData constructor')
        check_status(self.handle)

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def __add__(self, other):
        """Overloads + operator."""
        z = self.deep_copy()
        if isinstance(other, NiftiImageData):
            try_calling(pyreg.cReg_NiftiImageData_maths_im(z.handle, self.handle, other.handle, 0))
        else:
            try_calling(pyreg.cReg_NiftiImageData_maths_num(z.handle, self.handle, float(other), 0))
        check_status(z.handle)
        return z

    def __sub__(self, other):
        """Overloads - operator."""
        z = self.deep_copy()
        if isinstance(other, NiftiImageData):
            try_calling(pyreg.cReg_NiftiImageData_maths_im(z.handle, self.handle, other.handle, 1))
        else:
            try_calling(pyreg.cReg_NiftiImageData_maths_num(z.handle, self.handle, float(other), 1))
        check_status(z.handle)
        return z

    def __mul__(self, other):
        """Overloads * operator."""
        z = self.deep_copy()
        try_calling(pyreg.cReg_NiftiImageData_maths_num(z.handle, self.handle, float(other), 2))
        check_status(z.handle)
        return z

    def __eq__(self, other):
        """Overload comparison operator."""
        if not isinstance(other, NiftiImageData):
            raise AssertionError()
        if self.handle is None or other.handle is None:
            error('Cannot compare images as at least one is uninitialised')
        h = pyreg.cReg_NiftiImageData_equal(self.handle, other.handle)
        check_status(h, inspect.stack()[1])
        value = pyiutil.intDataFromHandle(h)
        pyiutil.deleteDataHandle(h)
        return value

    def __ne__(self, other):
        """Overload comparison operator."""
        return not self == other

    def write(self, filename, datatype=-1):
        """Save to file. See nifti1.h for datatypes (e.g., float (NIFTI_TYPE_FLOAT32) = 16).
        Image's original datatpye is used by default."""
        if self.handle is None:
            raise AssertionError()
        try_calling(pyreg.cReg_NiftiImageData_write(self.handle, filename, datatype))

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

    def get_dimensions(self):
        """Get dimensions. Returns nifti format.
        i.e., dim[0]=ndims, dim[1]=nx, dim[2]=ny,..."""
        if self.handle is None:
            raise AssertionError()
        dim = numpy.ndarray((8,), dtype=numpy.int32)
        try_calling(pyreg.cReg_NiftiImageData_get_dimensions(self.handle, dim.ctypes.data))
        return dim
    def get_voxel_sizes(self):
        """%Get voxel sizes. Returns nifti format.
        i.e., dim[0]=?, dim[1]=dx, dim[2]=dy,..."""
        if self.handle is None:
            raise AssertionError()
        out = numpy.ndarray((8,), dtype=numpy.float32)
        try_calling(pyreg.cReg_NiftiImageData_get_voxel_sizes(self.handle, out.ctypes.data))
        return out

    def fill(self, val):
        """Fill image with single value or numpy array."""
        if self.handle is None:
            raise AssertionError()
        if isinstance(val, numpy.ndarray):
            if val.dtype is numpy.dtype('float32'):
                v = val
            else:
                v = val.astype(numpy.float32)
            if not v.flags['F_CONTIGUOUS']:
                v = numpy.asfortranarray(v)
            try_calling(pyreg.cReg_NiftiImageData_fill_arr(self.handle, v.ctypes.data))
        elif isinstance(val, float):
            try_calling(pyreg.cReg_NiftiImageData_fill(self.handle, val))
        elif isinstance(val, int):
            try_calling(pyreg.cReg_NiftiImageData_fill(self.handle, float(val)))
        else:
            raise error('wrong fill value. Should be numpy.ndarray, float or int')

    def deep_copy(self):
        """Deep copy image."""
        if self.handle is None:
            raise AssertionError()
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
        try_calling(pyreg.cReg_NiftiImageData_deep_copy(image.handle, self.handle))
        return image

    def as_array(self):
        """Get data as numpy array."""
        if self.handle is None:
            raise AssertionError()
        dim = self.get_dimensions()
        dim = dim[1:dim[0]+1]
        array = numpy.ndarray(dim, dtype=numpy.float32, order='F')
        try_calling(pyreg.cReg_NiftiImageData_as_array(self.handle, array.ctypes.data))
        return numpy.ascontiguousarray(array)

    def get_original_datatype(self):
        """Get original image datatype (internally everything is converted to float)."""
        if self.handle is None:
            raise AssertionError()
        handle = pyreg.cReg_NiftiImageData_get_original_datatype(self.handle)
        check_status(handle)
        datatype = pyiutil.intDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return datatype

    def crop(self, min_, max_):
        """Crop image. Give minimum and maximum indices.
        Min and max indicies can be anywhere between (x,y,z) and (x,y,z,t,u,v,w).
        Use values of -1 for no change."""
        if len(min_) < 3 or len(min_) > 7:
            raise AssertionError("Min bounds should be at least (x,y,z), and up to (x,y,z,t,u,v,w)")
        if len(max_) < 3 or len(max_) > 7:
            raise AssertionError("Max bounds should be at least (x,y,z), and up to (x,y,z,t,u,v,w)")
        # Fill in any missing indices with -1's
        min_.extend([-1] * (7-len(min_)))
        max_.extend([-1] * (7-len(max_)))
        min_np = numpy.array(min_, dtype=numpy.int32)
        max_np = numpy.array(max_, dtype=numpy.int32)
        try_calling(pyreg.cReg_NiftiImageData_crop(self.handle, min_np.ctypes.data, max_np.ctypes.data))

    def print_header(self):
        """Print nifti header metadata."""
        vec = SIRF.DataHandleVector()
        vec.push_back(self.handle)
        try_calling(pyreg.cReg_NiftiImageData_print_headers(vec.handle))

    def same_object(self):
        """See DataContainer.same_object()."""
        return NiftiImageData()

    def set_voxel_spacing(self, spacing, interpolation_order):
        """Set the voxel spacing. Requires resampling image, and so interpolation order is required. 
        As per NiftyReg, interpolation_order can be either 0, 1 or 3 meaning nearest neighbor, linear or cubic spline interpolation."""
        if len(spacing) != 3:
            raise AssertionError("New spacing should be array of 3 numbers.")
            type(interpolation_order)
        try_calling(pyreg.cReg_NiftiImageData_set_voxel_spacing(self.handle, float(spacing[0]), float(spacing[1]), float(spacing[2]), int(interpolation_order)))

    def get_contains_nans(self):
        """Returns true if image contains any voxels with NaNs."""
        return parms.bool_par(self.handle, 'NiftiImageData', 'contains_nans')

    def normalise_zero_and_one(self):
        """Normalise image between 0 and 1."""
        try_calling(pyreg.cReg_NiftiImageData_normalise_zero_and_one(self.handle))
        check_status(self.handle)
    
    def standardise(self):
        """Standardise (subtract mean and divide by standard deviation)."""
        try_calling(pyreg.cReg_NiftiImageData_standardise(self.handle))
        check_status(self.handle)

    def get_inner_product(self, other):
        """Get inner product between two images. Must be same size."""
        if not isinstance(other, NiftiImageData):
            raise AssertionError()
        handle = pyreg.cReg_NiftiImageData_get_inner_product(self.handle, other.handle)
        check_status(handle)
        inner_product = pyiutil.floatDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return inner_product

    @staticmethod
    def print_headers(to_print):
        """Print nifti header metadata of one or multiple nifti images."""
        if not all(isinstance(n, NiftiImageData) for n in to_print):
            raise AssertionError()
        vec = SIRF.DataHandleVector()
        for n in to_print:
            vec.push_back(n.handle)
        try_calling(pyreg.cReg_NiftiImageData_print_headers(vec.handle))


class NiftiImageData3D(NiftiImageData):
    """
    Class for 3D nifti image data.
    """
    def __init__(self, src=None):
        self.handle = None
        self.name = 'NiftiImageData3D'
        if src is None:
            self.handle = pyreg.cReg_newObject(self.name)
        elif isinstance(src, str):
            self.handle = pyreg.cReg_objectFromFile(self.name, src)
        elif isinstance(src, SIRF.ImageData):
            # src is ImageData
            self.handle = pyreg.cReg_NiftiImageData3D_from_SIRFImageData(src.handle)
        else:
            raise error('Wrong source in NiftiImageData3D constructor')
        check_status(self.handle)

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)


class ImageData(NiftiImageData3D):
    """
    Alias class for nifti image data.
    """

class NiftiImageData3DTensor(NiftiImageData):
    """
    Class for 3D tensor nifti image data.
    """
    def __init__(self, src1=None, src2=None, src3=None):
        self.handle = None
        self.name = 'NiftiImageData3DTensor'
        if src1 is None:
            self.handle = pyreg.cReg_newObject(self.name)
        elif isinstance(src1, str):
            self.handle = pyreg.cReg_objectFromFile(self.name, src1)
        elif isinstance(src1, NiftiImageData3D) and isinstance(src2, NiftiImageData3D) and isinstance(src3, NiftiImageData3D):
            self.handle = pyreg.cReg_NiftiImageData3DTensor_construct_from_3_components(self.name, src1.handle,
                                                                                            src2.handle, src3.handle)
        else:
            raise error('Wrong source in NiftiImageData3DTensor constructor')
        check_status(self.handle)

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def write_split_xyz_components(self, filename, datatype=-1):
        """Save to file. See nifti1.h for datatypes (e.g., float (NIFTI_TYPE_FLOAT32) = 16).
        Image's original datatpye is used by default."""
        if self.handle is None:
            raise AssertionError()
        if not isinstance(filename, str):
            raise AssertionError()
        try_calling(pyreg.cReg_NiftiImageData3DTensor_write_split_xyz_components(self.handle, filename, datatype))

    def create_from_3D_image(self, src):
        """Create tensor/deformation/displacement field from 3D image."""
        if not isinstance(src, NiftiImageData3D):
            raise AssertionError()
        if src.handle is None:
            raise AssertionError()
        try_calling(pyreg.cReg_NiftiImageData3DTensor_create_from_3D_image(self.handle, src.handle))
        check_status(self.handle)

    def flip_component(self, dim):
        """Flip component of nu."""
        if 0 < dim or dim > 2:
            raise AssertionError("Dimension to flip should be between 0 and 2.")
        try_calling(pyreg.cReg_NiftiImageData3DTensor_flip_component(self.handle, dim))
        check_status(self.handle)


class NiftiImageData3DDisplacement(NiftiImageData3DTensor, _Transformation):
    """
    Class for 3D displacement nifti image data.
    """
    def __init__(self, src1=None, src2=None, src3=None):
        self.handle = None
        self.name = 'NiftiImageData3DDisplacement'
        if src1 is None:
            self.handle = pyreg.cReg_newObject(self.name)
        elif isinstance(src1, str):
            self.handle = pyreg.cReg_objectFromFile(self.name, src1)
        elif isinstance(src1, NiftiImageData3D) and isinstance(src2, NiftiImageData3D) and isinstance(src3, NiftiImageData3D):
            self.handle = pyreg.cReg_NiftiImageData3DTensor_construct_from_3_components(self.name, src1.handle,
                                                                                            src2.handle, src3.handle)
        elif isinstance(src1, NiftiImageData3DDeformation):
            self.handle = pyreg.cReg_NiftiImageData3DDisplacement_create_from_def(src1.handle)
        else:
            raise error('Wrong source in NiftiImageData3DDisplacement constructor')
        check_status(self.handle)

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)


class NiftiImageData3DDeformation(NiftiImageData3DTensor, _Transformation):
    """
    Class for 3D deformation nifti image data.
    """

    def __init__(self, src1=None, src2=None, src3=None):
        self.handle = None
        self.name = 'NiftiImageData3DDeformation'
        if src1 is None:
            self.handle = pyreg.cReg_newObject(self.name)
        elif isinstance(src1, str):
            self.handle = pyreg.cReg_objectFromFile(self.name, src1)
        elif isinstance(src1, NiftiImageData3D) and isinstance(src2, NiftiImageData3D) and isinstance(src3, NiftiImageData3D):
            self.handle = pyreg.cReg_NiftiImageData3DTensor_construct_from_3_components(self.name, src1.handle,
                                                                                            src2.handle, src3.handle)
        elif isinstance(src1, NiftiImageData3DDisplacement):
            self.handle = pyreg.cReg_NiftiImageData3DDeformation_create_from_disp(src1.handle)
        else:
            raise error('Wrong source in NiftiImageData3DDeformation constructor')
        check_status(self.handle)

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    @staticmethod
    def compose_single_deformation(trans, ref):
        """Compose transformations into single deformation."""
        if not isinstance(ref, NiftiImageData3D):
            raise AssertionError()
        if not all(isinstance(n, _Transformation) for n in trans):
            raise AssertionError()
        if len(trans) == 1:
            return trans[0].get_as_deformation_field(ref)
        # This is ugly. Store each type in a single string (need to do this because I can't get
        # virtual methods to work for multiple inheritance (deformation/displacement are both
        # nifti images and transformations).
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
        z.handle = pyreg.cReg_NiftiImageData3DDeformation_compose_single_deformation(
            ref.handle, types, vec.handle)
        check_status(z.handle)
        return z


class _Registration(ABC):
    """
    Abstract base class for registration.
    """
    def __init__(self):
        self.handle = None
        self.name = 'Registration'
        self.reference_image = None

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_parameter_file(self, filename):
        """Sets the parameter filename."""
        parms.set_char_par(self.handle, 'Registration', 'parameter_file', filename)

    def set_reference_image(self, reference_image):
        """Sets the reference image."""
        if not isinstance(reference_image, SIRF.ImageData):
            raise AssertionError()
        self.reference_image = reference_image
        parms.set_parameter(self.handle, 'Registration', 'reference_image', reference_image.handle)

    def set_floating_image(self, floating_image):
        """Sets the floating image."""
        if not isinstance(floating_image, SIRF.ImageData):
            raise AssertionError()
        parms.set_parameter(self.handle, 'Registration', 'floating_image', floating_image.handle)

    def set_reference_mask(self, reference_mask):
        """Sets the reference mask."""
        if not isinstance(reference_mask, SIRF.ImageData):
            raise AssertionError()
        parms.set_parameter(self.handle, 'Registration', 'reference_mask', reference_mask.handle)

    def set_floating_mask(self, floating_mask):
        """Sets the floating mask."""
        if not isinstance(floating_mask, SIRF.ImageData):
            raise AssertionError()
        parms.set_parameter(self.handle, 'Registration', 'floating_mask', floating_mask.handle)

    def get_output(self):
        """Gets the registered image."""
        output = self.reference_image.same_object()
        output.handle = parms.parameter_handle(self.handle, 'Registration', 'output')
        return output

    def process(self):
        """Run the registration"""
        try_calling(pyreg.cReg_Registration_process(self.handle))

    def get_deformation_field_forward(self):
        """Gets the forward deformation field image."""
        output = NiftiImageData3DDeformation()
        output.handle = pyreg.cReg_Registration_get_deformation_displacement_image(self.handle, 'forward_deformation')
        check_status(output.handle)
        return output

    def get_deformation_field_inverse(self):
        """Gets the inverse deformation field image."""
        output = NiftiImageData3DDeformation()
        output.handle = pyreg.cReg_Registration_get_deformation_displacement_image(self.handle, 'inverse_deformation')
        check_status(output.handle)
        return output

    def get_displacement_field_forward(self):
        """Gets the forward displacement field image."""
        output = NiftiImageData3DDisplacement()
        output.handle = pyreg.cReg_Registration_get_deformation_displacement_image(self.handle, 'forward_displacement')
        check_status(output.handle)
        return output

    def get_displacement_field_inverse(self):
        """Gets the inverse displacement field image."""
        output = NiftiImageData3DDisplacement()
        output.handle = pyreg.cReg_Registration_get_deformation_displacement_image(self.handle, 'inverse_displacement')
        check_status(output.handle)
        return output

    def set_parameter(self, par, arg1="", arg2=""):
        """Set string parameter. Check if any set methods match the method given by par.
        If so, set the value given by arg. Convert to float/int etc., as necessary.
        Up to 2 arguments, leave blank if unneeded. These are applied after parsing
        the parameter file."""
        try_calling(pyreg.cReg_Registration_set_parameter(self.handle, par, arg1, arg2))


class NiftyAladinSym(_Registration):
    """
    Registration using NiftyReg aladin.
    """
    def __init__(self):
        _Registration.__init__(self)
        self.name = 'NiftyAladinSym'
        self.handle = pyreg.cReg_newObject(self.name)
        check_status(self.handle)

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def get_transformation_matrix_forward(self):
        """Get forward transformation matrix."""
        if self.handle is None:
            raise AssertionError()
        tm = AffineTransformation()
        tm.handle = pyreg.cReg_NiftyAladin_get_TM(self.handle, 'forward')
        return tm

    def get_transformation_matrix_inverse(self):
        """Get inverse transformation matrix."""
        if self.handle is None:
            raise AssertionError()
        tm = AffineTransformation()
        tm.handle = pyreg.cReg_NiftyAladin_get_TM(self.handle, 'inverse')
        return tm

    @staticmethod
    def print_all_wrapped_methods():
        """Print all wrapped methods."""
        print("In C++, this class is templated. \"dataType\" corresponds to \"float\" for Matlab and python.")
        try_calling(pyreg.cReg_Registration_print_all_wrapped_methods('NiftyAladinSym'))


class NiftyF3dSym(_Registration):
    """
    Registration using NiftyReg f3d.
    """
    def __init__(self):
        _Registration.__init__(self)
        self.name = 'NiftyF3dSym'
        self.handle = pyreg.cReg_newObject(self.name)
        check_status(self.handle)

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

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
        print("In C++, this class is templated. \"dataType\" corresponds to \"float\" for Matlab and python.")
        try_calling(pyreg.cReg_Registration_print_all_wrapped_methods('NiftyF3dSym'))


class NiftyResample(object):
    """
    Resample using NiftyReg.
    """
    def __init__(self):
        self.name = 'NiftyResample'
        self.handle = pyreg.cReg_newObject(self.name)
        self.reference_image = None
        self.floating_image = None
        check_status(self.handle)

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_reference_image(self, reference_image):
        """Set reference image. This is the image that would be the reference if you were doing a forward transformation."""
        if not isinstance(reference_image, SIRF.ImageData):
            raise AssertionError()
        self.reference_image = reference_image
        parms.set_parameter(self.handle, self.name, 'reference_image', reference_image.handle)

    def set_floating_image(self, floating_image):
        """Set floating image. This is the image that would be the floating if you were doing a forward transformation."""
        if not isinstance(floating_image, SIRF.ImageData):
            raise AssertionError()
        self.floating_image = floating_image
        parms.set_parameter(self.handle, self.name, 'floating_image', floating_image.handle)

    def add_transformation(self, src):
        """Add transformation."""
        if isinstance(src, AffineTransformation):
            try_calling(pyreg.cReg_NiftyResample_add_transformation(self.handle, src.handle, 'affine'))
        elif isinstance(src, NiftiImageData3DDisplacement):
            try_calling(pyreg.cReg_NiftyResample_add_transformation(self.handle, src.handle, 'displacement'))
        elif isinstance(src, NiftiImageData3DDeformation):
            try_calling(pyreg.cReg_NiftyResample_add_transformation(self.handle, src.handle, 'deformation'))
        else:
            raise AssertionError()

    def set_interpolation_type(self, interp_type):
        """Set interpolation type. 0=nearest neighbour, 1=linear, 3=cubic, 4=sinc."""
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

    @deprecation.deprecated(deprecated_in="2.1.0",
                            details="Use forward(image) instead")
    def process(self):
        """Process."""
        try_calling(pyreg.cReg_NiftyResample_process(self.handle))

    def get_output(self):
        """Get output."""
        image = self.reference_image.same_object()
        image.handle = parms.parameter_handle(self.handle, self.name, 'output')
        check_status(image.handle)
        return image

    def forward(self, im1, im2=None):
        """Forward transformation.
        Usage:
            output = forward(input), OR
            forward(output,input)
        """
        # If usage was 'output = forward(input)'
        if im2 is None:
            output_im = self.reference_image.deep_copy()
            input_im = im1
        # If usage was 'forward(output, input)'
        else:
            output_im = im1
            input_im = im2
        # Check image validity
        if not isinstance(input_im, SIRF.ImageData):
            raise AssertionError()
        if not isinstance(output_im, SIRF.ImageData):
            raise AssertionError()
        # Forward
        try_calling(pyreg.cReg_NiftyResample_forward(output_im.handle, input_im.handle, self.handle))
        # If usage was 'output = forward(input)', we need to return the ouptut
        if im2 is None:
            return output_im

    def adjoint(self, im1, im2=None):
        """Adjoint transformation.
        Usage:
            output = adjoint(input), OR
            adjoint(output,input)
        """
        # If usage was 'output = adjoint(input)'
        if im2 is None:
            output_im = self.floating_image.deep_copy()
            input_im = im1
        # If usage was 'adjoint(output, input)'
        else:
            output_im = im1
            input_im = im2
        # Check image validity
        if not isinstance(input_im, SIRF.ImageData):
            raise AssertionError()
        if not isinstance(output_im, SIRF.ImageData):
            raise AssertionError()
        # Forward
        try_calling(pyreg.cReg_NiftyResample_adjoint(output_im.handle, input_im.handle, self.handle))
        # If usage was 'output = adjoint(input)', we need to return the ouptut
        if im2 is None:
            return output_im

    def backward(self, im1, im2=None):
        """Backward transformation. Alias of adjoint to align terms with AcquisitionModel's forward and backward."""
        return self.adjoint(im1, im2)


class ImageWeightedMean(object):
    """
    Class for performing weighted mean of images.
    """

    def __init__(self):
        self.name = 'ImageWeightedMean'
        self.handle = pyreg.cReg_newObject(self.name)
        check_status(self.handle)

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def add_image(self, image, weight):
        """Add an image (filename or NiftiImageData) and its corresponding weight."""
        if isinstance(image, NiftiImageData):
            try_calling(pyreg.cReg_ImageWeightedMean_add_image(self.handle, image.handle, weight))
        elif isinstance(image, str):
            try_calling(pyreg.cReg_ImageWeightedMean_add_image_filename(self.handle, image, weight))
        else:
            raise error("sirf.Reg.ImageWeightedMean.add_image: image must be NiftiImageData or filename.")

    def process(self):
        """Process."""
        try_calling(pyreg.cReg_ImageWeightedMean_process(self.handle))

    def get_output(self):
        """Get output."""
        image = NiftiImageData()
        image.handle = parms.parameter_handle(self.handle, self.name, 'output')
        check_status(image.handle)
        return image


class AffineTransformation(_Transformation):
    """
    Class for affine transformations.
    """
    def __init__(self, src=None, quat=None):
        self.handle = None
        self.name = 'AffineTransformation'
        if src is None:
            self.handle = pyreg.cReg_newObject(self.name)
        elif isinstance(src, str):
            self.handle = pyreg.cReg_objectFromFile(self.name, src)
        elif isinstance(src, numpy.ndarray) and quat is None:
            if src.shape != (4, 4):
                raise AssertionError()
            # Need to transpose relative to MATLAB
            trans = numpy.zeros((4, 4), dtype=numpy.float32)
            for i in range(4):
                for j in range(4):
                    trans[i,j] = src[j,i]
                self.handle = pyreg.cReg_AffineTransformation_construct_from_TM(trans.ctypes.data)
        elif isinstance(src, numpy.ndarray) and quat is not None and isinstance(quat, Quaternion):
            self.handle = pyreg.cReg_AffineTransformation_construct_from_trans_and_quaternion(src.ctypes.data, quat.handle)
        else:
            raise error('AffineTransformation accepts no args, filename, 4x4 array or translation with quaternion.')
        check_status(self.handle)

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def __eq__(self, other):
        """Overload comparison operator."""
        if not isinstance(other, AffineTransformation):
            raise AssertionError()
        h = pyreg.cReg_AffineTransformation_equal(self.handle, other.handle)
        check_status(h, inspect.stack()[1])
        value = pyiutil.intDataFromHandle(h)
        pyiutil.deleteDataHandle(h)
        return value

    def __ne__(self, other):
        """Overload comparison operator."""
        return not self == other

    def __mul__(self, other):
        """Overload multiplication operator."""
        if not isinstance(other, AffineTransformation):
            raise AssertionError()
        mat = AffineTransformation()
        mat.handle = pyreg.cReg_AffineTransformation_mul(self.handle, other.handle)
        check_status(mat.handle)
        return mat

    def deep_copy(self):
        """Deep copy."""
        if self.handle is None:
            raise AssertionError()
        mat = AffineTransformation()
        mat.handle = pyreg.cReg_AffineTransformation_deep_copy(self.handle)
        check_status(mat.handle)
        return mat

    def write(self, filename):
        """Save to file."""
        if self.handle is None:
            raise AssertionError()
        try_calling(pyreg.cReg_AffineTransformation_write(self.handle, filename))

    def get_determinant(self):
        """Get determinant."""
        return parms.float_par(self.handle, self.name, 'determinant')

    def as_array(self):
        """Get forward transformation matrix."""
        if self.handle is None:
            raise AssertionError()
        tm = numpy.ndarray((4, 4), dtype=numpy.float32)
        try_calling(pyreg.cReg_AffineTransformation_as_array(self.handle, tm.ctypes.data))
        return tm

    def get_inverse(self):
        """Get inverse matrix."""
        tm = AffineTransformation()
        tm.handle = pyreg.cReg_AffineTransformation_get_inverse(self.handle)
        check_status(tm.handle)
        return tm

    def get_Euler_angles(self):
        """Get Euler angles of transformation matrix (XYZ)."""
        if self.handle is None:
            raise AssertionError()
        eul = numpy.ndarray(3, dtype=numpy.float32)
        try_calling(pyreg.cReg_AffineTransformation_get_Euler_angles(self.handle, eul.ctypes.data))
        return eul

    def get_quaternion(self):
        """Get quaternion."""
        if self.handle is None:
            raise AssertionError()
        quat_zeros = numpy.array([0., 0., 0., 0.],dtype=numpy.float32)
        quat = Quaternion(quat_zeros)
        quat.handle = pyreg.cReg_AffineTransformation_get_quaternion(self.handle)
        check_status(quat.handle)
        return quat

    @staticmethod
    def get_identity():
        """Get identity matrix."""
        mat = AffineTransformation()
        mat.handle = pyreg.cReg_AffineTransformation_get_identity()
        return mat

    @staticmethod
    def get_average(to_average):
        """Get average of transformations."""
        if not all(isinstance(n, AffineTransformation) for n in to_average):
            raise AssertionError("AffineTransformation:get_average() input list should only contain AffineTransformations.")
        tm = AffineTransformation()
        vec = SIRF.DataHandleVector()
        for n in to_average:
            vec.push_back(n.handle)
        tm.handle = pyreg.cReg_AffineTransformation_get_average(vec.handle)
        check_status(tm.handle)
        return tm


class Quaternion(object):
    """
    Class for quaternions.
    """
    def __init__(self, src=None):
        self.handle = None
        self.name = 'Quaternion'
        if isinstance(src, numpy.ndarray):
            if src.size != 4:
                raise AssertionError('Quaternion constructor from numpy array is wrong size.')
            self.handle = pyreg.cReg_Quaternion_construct_from_array(src.ctypes.data)
        elif isinstance(src, AffineTransformation):
            self.handle = pyreg.cReg_Quaternion_construct_from_AffineTransformation(src.handle)
        else:
            raise error('Wrong source in quaternion constructor')
        check_status(self.handle)

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def as_array(self):
        """Get quaternion as array."""
        if self.handle is None:
            raise AssertionError()
        arr = numpy.ndarray(4, dtype=numpy.float32)
        try_calling(pyreg.cReg_Quaternion_as_array(self.handle, arr.ctypes.data))
        return arr

    @staticmethod
    def get_average(to_average):
        """Get average of quaternions."""
        if not all(isinstance(n, Quaternion) for n in to_average):
            raise AssertionError()
        quat_zeros = numpy.array([0., 0., 0., 0.],dtype=numpy.float32)
        quat = Quaternion(quat_zeros)
        if not all(isinstance(n, Quaternion) for n in to_average):
            raise AssertionError()
        vec = SIRF.DataHandleVector()
        for n in to_average:
            vec.push_back(n.handle)
        quat.handle = pyreg.cReg_Quaternion_get_average(vec.handle)
        return quat
