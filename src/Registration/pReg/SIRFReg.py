""" 
Object-Oriented wrap for the cSIRFReg-to-Python interface pysirfreg.py
"""

# CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
# Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
# Copyright 2015 - 2018 University College London
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
import pysirfreg
import pSTIR
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


def _setParameter(hs, set, par, hv, stack=None):
    # try_calling(pysirfreg.cSIRFReg_setParameter(hs, set, par, hv))
    if stack is None:
        stack = inspect.stack()[1]
    h = pysirfreg.cSIRFReg_setParameter(hs, set, par, hv)
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
    h = pysirfreg.cSIRFReg_parameter(handle, set, par)
    check_status(h, inspect.stack()[1])
    value = pyiutil.charDataFromHandle(h)
    pyiutil.deleteDataHandle(h)
    return value


def _int_par(handle, set, par):
    h = pysirfreg.cSIRFReg_parameter(handle, set, par)
    check_status(h, inspect.stack()[1])
    value = pyiutil.intDataFromHandle(h)
    pyiutil.deleteDataHandle(h)
    return value


def _float_par(handle, set, par):
    h = pysirfreg.cSIRFReg_parameter(handle, set, par)
    check_status(h, inspect.stack()[1])
    value = pyiutil.floatDataFromHandle(h)
    pyiutil.deleteDataHandle(h)
    return value


def _float_pars(handle, set, par, n):
    h = pysirfreg.cSIRFReg_parameter(handle, set, par)
    check_status(h)
    value = ()
    for i in range(n):
        value += (pyiutil.floatDataItemFromHandle(h, i), )
    pyiutil.deleteDataHandle(h)
    return value


def _getParameterHandle(hs, set, par):
    handle = pysirfreg.cSIRFReg_parameter(hs, set, par)
    check_status(handle, inspect.stack()[1])
    return handle


def _tmp_filename():
    return repr(int(1000*time.time()))
###########################################################


class MessageRedirector:
    """
    Class for SIRFReg printing redirection to files/stdout/stderr.
    """
    def __init__(self, info=None, warn='stdout', errr='stdout'):
        """
        Creates MessageRedirector object that redirects SIRFReg's ouput
        produced by info(), warning() and error(0 functions to destinations
        specified respectively by info, warn and err arguments.
        The argument values other than None, stdout, stderr, cout and cerr
        are interpreted as filenames.
        None and empty string value suppresses printing.
        """
        if info is None:
            info = ''
        if type(info) is not type(' '):
            raise error('wrong info argument for MessageRedirector constructor')
        elif info in {'stdout', 'stderr', 'cout', 'cerr'}:
            self.info = pysirfreg.newTextPrinter(info)
            self.info_case = 0
        else:
            self.info = pysirfreg.newTextWriter(info)
            self.info_case = 1
        pysirfreg.openChannel(0, self.info)

        if warn is None:
            warn = ''
        if type(warn) is not type(' '):
            raise error('wrong warn argument for MessageRedirector constructor')
        elif warn in {'stdout', 'stderr', 'cout', 'cerr'}:
            self.warn = pysirfreg.newTextPrinter(warn)
            self.warn_case = 0
        else:
            self.warn = pysirfreg.newTextWriter(warn)
            self.warn_case = 1
        pysirfreg.openChannel(1, self.warn)

        if errr is None:
            errr = ''
        if type(errr) is not type(' '):
            raise error('wrong errr argument for MessageRedirector constructor')
        elif errr in {'stdout', 'stderr', 'cout', 'cerr'}:
            self.errr = pysirfreg.newTextPrinter(errr)
            self.errr_case = 0
        else:
            self.errr = pysirfreg.newTextWriter(errr)
            self.errr_case = 1
        pysirfreg.openChannel(2, self.errr)

    def __del__(self):
        if self.info_case == 0:
            try_calling(pysirfreg.deleteTextPrinter(self.info))
        else:
            try_calling(pysirfreg.deleteTextWriter(self.info))
        pysirfreg.closeChannel(0, self.info)
        if self.warn_case == 0:
            try_calling(pysirfreg.deleteTextPrinter(self.warn))
        else:
            try_calling(pysirfreg.deleteTextWriter(self.warn))
        pysirfreg.closeChannel(1, self.warn)
        if self.errr_case == 0:
            try_calling(pysirfreg.deleteTextPrinter(self.errr))
        else:
            try_calling(pysirfreg.deleteTextWriter(self.errr))
        pysirfreg.closeChannel(2, self.errr)
###########################################################


def do_nifti_images_match(im1, im2, accuracy_percentage_of_max):
    """Do nifti images match?"""
    assert isinstance(im1, ImageData)
    assert isinstance(im2, ImageData)
    handle = pysirfreg.cSIRFReg_do_nifti_images_match(im1.handle, im2.handle, accuracy_percentage_of_max)
    check_status(handle)
    r = pyiutil.intDataFromHandle(handle)
    pyiutil.deleteDataHandle(handle)
    return r


def dump_nifti_info(to_dump):
    """Dump metadata of one or multiple (up to 5) nifti images."""
    if isinstance(to_dump, str):
        try_calling(pysirfreg.cSIRFReg_dump_nifti_info_filename(to_dump))
    elif isinstance(to_dump, ImageData):
        try_calling(pysirfreg.cSIRFReg_dump_nifti_info_im1(to_dump.handle))
    elif all(isinstance(n, ImageData) for n in to_dump):
        if len(to_dump) == 1:
            try_calling(pysirfreg.cSIRFReg_dump_nifti_info_im1(to_dump[0].handle))
        elif len(to_dump) == 2:
            try_calling(pysirfreg.cSIRFReg_dump_nifti_info_im2(to_dump[0].handle, to_dump[1].handle))
        elif len(to_dump) == 3:
            try_calling(pysirfreg.cSIRFReg_dump_nifti_info_im3(to_dump[0].handle, to_dump[1].handle, to_dump[2].handle))
        elif len(to_dump) == 4:
            try_calling(pysirfreg.cSIRFReg_dump_nifti_info_im4(to_dump[0].handle, to_dump[1].handle, to_dump[2].handle,
                                                               to_dump[3].handle))
        elif len(to_dump) == 5:
            try_calling(pysirfreg.cSIRFReg_dump_nifti_info_im5(to_dump[0].handle, to_dump[1].handle, to_dump[2].handle,
                                                               to_dump[3].handle, to_dump[4].handle))
        else:
            raise error('dump_nifti_info only implemented for up to 5 images.')

    else:
        raise error('dump_nifti_info requires filename, SIRFImageData or a list of SIRFImageData.')


def get_matrix(filename=None):
    """Get 4x4 matrix. If str is given, read from file. Else, return identity."""
    if filename is None:
        return numpy.identity(4, dtype=numpy.float32)
    assert isinstance(filename, str)
    tm = numpy.ndarray((4, 4), dtype=numpy.float32)
    try_calling(pysirfreg.cSIRFReg_SIRFReg_open_TM(filename, tm.ctypes.data))
    return tm


def compose_transformations_into_single_deformation(trans, ref):
    """Compose up to transformations into single deformation."""
    assert isinstance(ref, ImageData)
    assert all(isinstance(n, _Transformation) for n in trans)
    if len(trans) == 1:
        return trans[0]
    z = TransformationDeformation()
    if len(trans) == 2:
        z.handle = pysirfreg.cSIRFReg_compose_transformations_into_single_deformation2(
            ref.handle, trans[0].handle, trans[1].handle)
    elif len(trans) == 3:
        z.handle = pysirfreg.cSIRFReg_compose_transformations_into_single_deformation3(
            ref.handle, trans[0].handle, trans[1].handle, trans[2].handle)
    elif len(trans) == 4:
        z.handle = pysirfreg.cSIRFReg_compose_transformations_into_single_deformation4(
            ref.handle, trans[0].handle, trans[1].handle, trans[2].handle, trans[3].handle)
    elif len(trans) == 5:
        z.handle = pysirfreg.cSIRFReg_compose_transformations_into_single_deformation5(
            ref.handle, trans[0].handle, trans[1].handle, trans[2].handle, trans[3].handle, trans[4].handle)
    else:
        raise error('compose_transformations_into_single_deformation only implemented for up to 5 transformations.')
    check_status(z.handle)
    return z


class ImageData:
    """
    Class for image data.
    """
    def __init__(self, src=None):
        self.handle = None
        self.name = 'SIRFImageData'
        if src is None:
            self.handle = pysirfreg.cSIRFReg_newObject(self.name)
        elif isinstance(src, str):
            self.handle = pysirfreg.cSIRFReg_objectFromFile(self.name, src)
        elif isinstance(src, pSTIR.ImageData):
            # src is ImageData
            self.handle = pysirfreg.cSIRFReg_SIRFImageData_from_PETImageData(src.handle)
        else:
            raise error('Wrong source in ImageData constructor')
        check_status(self.handle)

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def __add__(self, other):
        """
        Overloads + for image data.
        """
        assert_validities(self, other)
        z = ImageData()
        z.handle = pysirfreg.cSIRFReg_SIRFImageData_maths(self.handle, other.handle, 1)
        check_status(z.handle)
        return z

    def __sub__(self, other):
        """
        Overloads - for image data.
        """
        assert_validities(self, other)
        z = ImageData()
        z.handle = pysirfreg.cSIRFReg_SIRFImageData_maths(self.handle, other.handle, -1)
        check_status(z.handle)
        return z

    def save_to_file(self, filename):
        """Save to file."""
        assert self.handle is not None
        try_calling(pysirfreg.cSIRFReg_SIRFImageData_save_to_file(self.handle, filename))

    def get_max(self):
        """Get max."""
        return _float_par(self.handle, 'SIRFImageData', 'max')

    def get_min(self):
        """Get min."""
        return _float_par(self.handle, 'SIRFImageData', 'min')

    def get_sum(self):
        """Get sum."""
        return _float_par(self.handle, 'SIRFImageData', 'sum')

    def get_dimensions(self):
        """Get dimensions. Returns nifti format.
        i.e., dim[0]=ndims, dim[1]=nx, dim[2]=ny,..."""
        assert self.handle is not None
        dim = numpy.ndarray((8,), dtype=numpy.int32)
        try_calling(pysirfreg.cSIRFReg_SIRFImageData_get_dimensions(self.handle, dim.ctypes.data))
        return dim

    def copy_data_to(self, pet_image):
        """Fill the STIRImageData with the values from SIRFImageData."""
        assert self.handle is not None
        assert isinstance(pet_image, pSTIR.ImageData)
        try_calling(pysirfreg.cSIRFReg_SIRFImageData_copy_data_to(self.handle, pet_image.handle))

    def fill(self, val):
        """Fill image with single value."""
        assert self.handle is not None
        try_calling(pysirfreg.cSIRFReg_SIRFImageData_fill(self.handle, val))

    def deep_copy(self):
        """Deep copy image."""
        assert self.handle is not None
        image = ImageData()
        image.handle = pysirfreg.cSIRFReg_SIRFImageData_deep_copy(self.handle)
        return image

    def as_array(self):
        """Get data as numpy array."""
        assert self.handle is not None
        dim = self.get_dimensions()
        dim = dim[1:dim[0]+1]
        array = numpy.ndarray(dim, dtype=numpy.float32)
        try_calling(pysirfreg.cSIRFReg_SIRFImageData_get_data(self.handle, array.ctypes.data))
        return array


class ImageDataDeformation(ImageData):
    """
    Class for deformation/displacement image data.
    """
    def __init__(self, src=None):
        self.handle = None
        self.name = 'SIRFImageDataDeformation'
        if src is None:
            self.handle = pysirfreg.cSIRFReg_newObject(self.name)
        elif isinstance(src, str):
            self.handle = pysirfreg.cSIRFReg_objectFromFile(self.name, src)
        else:
            raise error('Wrong source in ImageDataDeformation constructor')
        check_status(self.handle)

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def save_to_file_split_xyz_components(self, filename):
        """Save to file."""
        assert self.handle is not None
        assert isinstance(filename, str)
        try_calling(pysirfreg.cSIRFReg_SIRFImageDataDeformation_save_to_file_split_xyz_components(self.handle, filename))

    def create_from_3D_image(self, src):
        """Create deformation/displacement field from 3D image."""
        assert self.handle is not None
        assert src.handle is not None
        assert isinstance(src, ImageData)
        try_calling(pysirfreg.cSIRFReg_SIRFImageDataDeformation_create_from_3D_image(self.handle, src.handle))
        check_status(self.handle)

    def deep_copy(self):
        """Deep copy image."""
        assert self.handle is not None
        image = ImageDataDeformation()
        image.handle = pysirfreg.cSIRFReg_SIRFImageDataDeformation_deep_copy(self.handle)
        return image


class _SIRFReg(ABC):
    """
    Abstract base class for registration.
    """
    def __init__(self):
        self.handle = None
        self.name = 'SIRFReg'

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_parameter_file(self, filename):
        """Sets the parameter filename."""
        _set_char_par(self.handle, 'SIRFReg', 'parameter_file', filename)

    def set_reference_image(self, src):
        """Sets the reference image."""
        assert isinstance(src, ImageData)
        _setParameter(self.handle, 'SIRFReg', 'reference_image', src.handle)

    def set_floating_image(self, src):
        """Sets the floating image."""
        assert isinstance(src, ImageData)
        _setParameter(self.handle, 'SIRFReg', 'floating_image', src.handle)

    def get_output(self):
        """Gets the registered image."""
        output = ImageData()
        output.handle = pysirfreg.cSIRFReg_parameter(self.handle, 'SIRFReg', 'output')
        check_status(output.handle)
        return output


    def update(self):
        """Run the registration"""
        try_calling(pysirfreg.cSIRFReg_SIRFReg_update(self.handle))

    def get_deformation_field_fwrd(self):
        """Gets the forward deformation field image."""
        output = ImageDataDeformation()
        output.handle = pysirfreg.cSIRFReg_SIRFReg_get_deformation_displacement_image(self.handle, 'fwrd_deformation')
        check_status(output.handle)
        return output

    def get_deformation_field_back(self):
        """Gets the backwards deformation field image."""
        output = ImageDataDeformation()
        output.handle = pysirfreg.cSIRFReg_SIRFReg_get_deformation_displacement_image(self.handle, 'back_deformation')
        check_status(output.handle)
        return output

    def get_displacement_field_fwrd(self):
        """Gets the forward displacement field image."""
        output = ImageDataDeformation()
        output.handle = pysirfreg.cSIRFReg_SIRFReg_get_deformation_displacement_image(self.handle, 'fwrd_displacement')
        check_status(output.handle)
        return output

    def get_displacement_field_back(self):
        """Gets the backwards displacement field image."""
        output = ImageDataDeformation()
        output.handle = pysirfreg.cSIRFReg_SIRFReg_get_deformation_displacement_image(self.handle, 'back_displacement')
        check_status(output.handle)
        return output


class NiftyAladinSym(_SIRFReg):
    """
    Registration using NiftyReg aladin.
    """
    def __init__(self):
        self.name = 'SIRFRegNiftyAladinSym'
        self.handle = pysirfreg.cSIRFReg_newObject(self.name)
        check_status(self.handle)

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def save_transformation_matrix_fwrd(self, filename):
        """Save forward transformation matrix."""
        assert isinstance(filename, str)
        try_calling(pysirfreg.cSIRFReg_SIRFRegNiftyAladinSym_save_transformation_matrix(self.handle, filename, 'fwrd'))

    def save_transformation_matrix_back(self, filename):
        """Save backward transformation matrix."""
        assert isinstance(filename, str)
        try_calling(pysirfreg.cSIRFReg_SIRFRegNiftyAladinSym_save_transformation_matrix(self.handle, filename, 'back'))

    def get_transformation_matrix_fwrd(self):
        """Get forward transformation matrix."""
        assert self.handle is not None
        tm = numpy.ndarray((4, 4), dtype=numpy.float32)
        try_calling(pysirfreg.cSIRFReg_SIRFReg_get_TM(self.handle, tm.ctypes.data, 'fwrd'))
        return tm

    def get_transformation_matrix_back(self):
        """Get backwards transformation matrix."""
        assert self.handle is not None
        tm = numpy.ndarray((4, 4), dtype=numpy.float32)
        try_calling(pysirfreg.cSIRFReg_SIRFReg_get_TM(self.handle, tm.ctypes.data, 'back'))
        return tm


class NiftyF3dSym(_SIRFReg):
    """
    Registration using NiftyReg f3d.
    """
    def __init__(self):
        self.name = 'SIRFRegNiftyF3dSym'
        self.handle = pysirfreg.cSIRFReg_newObject(self.name)
        check_status(self.handle)

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_floating_time_point(self, floating_time_point):
        """Set floating time point."""
        _set_int_par(self.handle, self.name, 'floating_time_point', floating_time_point)

    def set_reference_time_point(self, reference_time_point):
        """Set reference time point."""
        _set_int_par(self.handle, self.name, 'reference_time_point', reference_time_point)

    def set_initial_affine_transformation(self, initial_affine_transformation):
        """Set initial affine transformation."""
        _set_char_par(self.handle, self.name, 'initial_affine_transformation', initial_affine_transformation)


class NiftyResample:
    """
    Resample using NiftyReg.
    """
    def __init__(self):
        self.name = 'SIRFRegNiftyResample'
        self.handle = pysirfreg.cSIRFReg_newObject(self.name)
        check_status(self.handle)

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_reference_image(self, reference_image):
        """Set reference image."""
        assert isinstance(reference_image, ImageData)
        _setParameter(self.handle, self.name, 'reference_image', reference_image.handle)

    def set_floating_image(self, floating_image):
        """Set floating image."""
        assert isinstance(floating_image, ImageData)
        _setParameter(self.handle, self.name, 'floating_image', floating_image.handle)

    def add_transformation_affine(self, src):
        """Add affine transformation."""
        assert isinstance(src, TransformationAffine)
        try_calling(pysirfreg.cSIRFReg_SIRFRegNiftyResample_add_transformation(self.handle, src.handle, 'affine'))

    def add_transformation_disp(self, src):
        """Add displacement field."""
        sys.stderr.write('\n\n\nim here1\n\n\n\n')
        assert isinstance(src, TransformationDisplacement)
        sys.stderr.write('\n\n\nim here2\n\n\n\n')
        try_calling(pysirfreg.cSIRFReg_SIRFRegNiftyResample_add_transformation(self.handle, src.handle, 'displacement'))
        sys.stderr.write('\n\n\nim here3\n\n\n\n')

    def add_transformation_def(self, src):
        """Add deformation field."""
        assert isinstance(src, TransformationDeformation)
        try_calling(pysirfreg.cSIRFReg_SIRFRegNiftyResample_add_transformation(self.handle, src.handle, 'deformation'))

    def set_interpolation_type(self, type):
        """Set interpolation type. 0=nearest neighbour, 1=linear, 3=cubic, 4=sinc."""
        assert isinstance(type, int)
        _set_int_par(self.handle, self.name, 'interpolation_type', type)

    def set_interpolation_type_to_nearest_neighbour(self):
        """Set interpolation type to nearest neighbour."""
        _set_int_par(self.handle, self.name, 'interpolation_type', 0)

    def set_interpolation_type_to_linear(self):
        """Set interpolation type to linear."""
        _set_int_par(self.handle, self.name, 'interpolation_type', 1)

    def set_interpolation_type_to_cubic_spline(self):
        """Set interpolation type to cubic spline."""
        _set_int_par(self.handle, self.name, 'interpolation_type', 3)

    def set_interpolation_type_to_sinc(self):
        """Set interpolation type to sinc."""
        _set_int_par(self.handle, self.name, 'interpolation_type', 4)

    def update(self):
        """Update."""
        try_calling(pysirfreg.cSIRFReg_SIRFRegNiftyResample_update(self.handle))

    def get_output(self):
        """Get output."""
        image = ImageData()
        image.handle = _getParameterHandle(self.handle, self.name, 'output')
        check_status(image.handle)
        return image


class ImageWeightedMean3D:
    """
    Class for performing weighted mean of images.
    """

    def __init__(self):
        self.name = 'SIRFRegImageWeightedMean3D'
        self.handle = pysirfreg.cSIRFReg_newObject(self.name)
        check_status(self.handle)

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def add_image(self, image, weight):
        """Add an image (filename or SIRFImageData) and its corresponding weight."""
        if isinstance(image, ImageData):
            try_calling(pysirfreg.cSIRFReg_SIRFRegImageWeightedMean3D_add_image(self.handle, image.handle, weight))
        elif isinstance(image, str):
            try_calling(pysirfreg.cSIRFReg_SIRFRegImageWeightedMean3D_add_image_filename(self.handle, image, weight))
        else:
            raise error("pSIRFReg.ImageWeightedMean3D.add_image: image must be SIRFImageData or filename.")

    def update(self):
        """Update."""
        try_calling(pysirfreg.cSIRFReg_SIRFRegImageWeightedMean3D_update(self.handle))

    def get_output(self):
        """Get output."""
        image = ImageData()
        image.handle = _getParameterHandle(self.handle, self.name, 'output')
        check_status(image.handle)
        return image


class ImageWeightedMean4D:
    """
    Class for performing weighted mean of deformation/displacement field images.
    """
    def __init__(self):
        self.name = 'SIRFRegImageWeightedMean4D'
        self.handle = pysirfreg.cSIRFReg_newObject(self.name)
        check_status(self.handle)

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def add_image(self, image, weight):
        """Add an image (filename or SIRFImageDataDeformation) and its corresponding weight."""
        if isinstance(image, ImageData):
            try_calling(pysirfreg.cSIRFReg_SIRFRegImageWeightedMean4D_add_image(self.handle, image.handle, weight))
        elif isinstance(image, str):
            try_calling(pysirfreg.cSIRFReg_SIRFRegImageWeightedMean4D_add_image_filename(self.handle, image, weight))
        else:
            raise error("pSIRFReg.ImageWeightedMean4D.add_image: image must be SIRFImageDataDeformation or filename.")

    def update(self):
        """Update."""
        try_calling(pysirfreg.cSIRFReg_SIRFRegImageWeightedMean4D_update(self.handle))

    def get_output(self):
        """Get output."""
        image = ImageData()
        image.handle = _getParameterHandle(self.handle, self.name, 'output')
        check_status(image.handle)
        return image


class _Transformation(ABC):
    """
    Abstract base class for transformations.
    """
    def __init__(self):
        self.handle = None
        self.name = 'SIRFRegTransformation'

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def get_as_deformation_field(self, ref):
        """Get any type of transformation as a deformation field.
        This is useful for joining them together. Require a reference
        image for converting transformation matrices to deformations."""
        assert isinstance(ref, ImageData)
        output = ImageDataDeformation()
        output.handle = pysirfreg.cSIRFReg_SIRFRegTransformation_get_as_deformation_field(self.handle, ref.handle)
        check_status(output.handle)
        return output


class TransformationAffine(_Transformation):
    """
    Class for affine transformations.
    """
    def __init__(self, src=None):
        self.name = 'SIRFRegTransformationAffine'
        if src is None:
            self.handle = pysirfreg.cSIRFReg_newObject(self.name)
        elif isinstance(src, str):
            self.handle = pysirfreg.cSIRFReg_objectFromFile(self.name, src)
        elif isinstance(src, numpy.ndarray):
            assert(src.shape == (4, 4))
            self.handle = pysirfreg.cSIRFReg_SIRFRegTransformationAffine_construct_from_TM(src.ctypes.data)
        else:
            raise error('Wrong source in affine transformation constructor')
        check_status(self.handle)

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)


class TransformationDisplacement(_Transformation):
    """
    Class for displacement transformations.
    """
    def __init__(self, src=None):
        self.name = 'SIRFRegTransformationDisplacement'
        if src is None:
            self.handle = pysirfreg.cSIRFReg_newObject(self.name)
        elif isinstance(src, str):
            self.handle = pysirfreg.cSIRFReg_objectFromFile(self.name, src)
        elif isinstance(src, ImageDataDeformation):
            self.handle = pysirfreg.cSIRFReg_SIRFRegTransformationDisplacement_construct_from_SIRFImageDataDeformation(
                src.handle)
        else:
            raise error('Wrong source in displacement transformation constructor')
        check_status(self.handle)

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)


class TransformationDeformation(_Transformation):
    """
    Class for deformation transformations.
    """
    def __init__(self, src=None):
        self.name = 'SIRFRegTransformationDeformation'
        if src is None:
            self.handle = pysirfreg.cSIRFReg_newObject(self.name)
        elif isinstance(src, str):
            self.handle = pysirfreg.cSIRFReg_objectFromFile(self.name, src)
        elif isinstance(src, ImageDataDeformation):
            self.handle = pysirfreg.cSIRFReg_SIRFRegTransformationDeformation_construct_from_SIRFImageDataDeformation(
                src.handle)
        else:
            raise error('Wrong source in deformation transformation constructor')
        check_status(self.handle)

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

