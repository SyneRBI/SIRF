""" 
Object-Oriented wrap for the cSIRFReg-to-Python interface pysirfreg.py
"""

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
## Copyright 2015 - 2018 University College London
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


def _setParameter(hs, set, par, hv, stack = None):
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


class ImageData:
    """
    Class for image data.
    """
    def __init__(self, src = None):
        self.handle = None
        self.name = 'SIRFImageData'
        if src is None:
            self.handle = pysirfreg.cSIRFReg_newObject(self.name)
        elif isinstance(src, str):
            self.handle = pysirfreg.cSIRFReg_objectFromFile(self.name, src)
        elif isinstance(src, pSTIR.ImageData):
            # src is ImageData
            self.handle = pysirfreg.cSIRFReg_SIRFImageData_from_PETImageData\
                (src.handle)
        else:
            raise error('Wrong source in ImageData constructor')
        check_status(self.handle)

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def save_to_file(self, filename):
        """Save to file."""
        assert self.handle is not None
        try_calling(pysirfreg.cSIRFReg_SIRFImageData_save_to_file(self.handle, filename))

    def get_max(self):
        """Get max."""
        return _float_par(self.handle, self.name, 'max')

    def get_min(self):
        """Get min."""
        return _float_par(self.handle, self.name, 'min')

    def copy_data_to(self, pet_image):
        """Fill the STIRImageData with the values from SIRFImageData."""
        assert self.handle is not None
        assert isinstance(pet_image, pSTIR.ImageData)
        try_calling(pysirfreg.cSIRFReg_SIRFImageData_copy_data_to(self.handle, pet_image.handle))

    def fill(self, val):
        """Fill image with single value."""
        assert self.handle is not None
        try_calling(pysirfreg.cSIRFReg_SIRFImageData_fill(self.handle, val))


class ImageDataDeformation(ImageData):
    """
    Class for deformation/displacement image data.
    """
    def __init__(self, src = None):
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

    def save_to_file(self, filename, split_xyz):
        """Save to file."""
        assert self.handle is not None
        assert isinstance(filename, str)
        assert isinstance(split_xyz, bool)
        try_calling(pysirfreg.cSIRFReg_SIRFImageDataDeformation_save_to_file(self.handle, filename, split_xyz))

    def create_from_3D_image(self, src):
        """Create deformation/displacement field from 3D image."""
        assert self.handle is not None
        assert src.handle is not None
        assert isinstance(src, ImageData)
        try_calling(pysirfreg.cSIRFReg_SIRFImageDataDeformation_create_from_3D_image(self.handle, src.handle))


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

    def set_reference_image(self, input):
        """Sets the reference image."""
        assert isinstance(input, ImageData)
        _setParameter(self.handle, 'SIRFReg', 'reference_image', input.handle)

    def set_floating_image(self, input):
        """Sets the floating image."""
        assert isinstance(input, ImageData)
        _setParameter(self.handle, 'SIRFReg', 'floating_image', input.handle)

    def get_output(self):
        """Gets the registered image."""
        output = SIRFImageData()
        output.handle = pysirfreg.cSIRFReg_parameter(self.handle, 'SIRFReg', 'output')
        check_status(output.handle)
        return output

    def save_warped_image(self, filename):
        """Save warped image."""
        assert isinstance(filename, str)
        try_calling(pysirfreg.cSIRFReg_SIRFReg_save_image(self.handle, filename))

    def save_deformation_field_fwrd_image(self, filename, split_xyz):
        """Save forward deformation field image to file."""
        assert isinstance(filename, str)
        assert isinstance(split_xyz, bool)
        try_calling(pysirfreg.cSIRFReg_SIRFReg_save_deformation_displacement_image(self.handle, filename, 'fwrd_deformation', split_xyz))

    def save_deformation_field_back_image(self, filename, split_xyz):
        """Save backward deformation field image to file."""
        assert isinstance(filename, str)
        assert isinstance(split_xyz, bool)
        try_calling(pysirfreg.cSIRFReg_SIRFReg_save_deformation_displacement_image(self.handle, filename, 'back_deformation', split_xyz))

    def save_displacement_field_fwrd_image(self, filename, split_xyz):
        """Save forward displacement field image to file."""
        assert isinstance(filename, str)
        assert isinstance(split_xyz, bool)
        try_calling(pysirfreg.cSIRFReg_SIRFReg_save_deformation_displacement_image(self.handle, filename, 'fwrd_displacement', split_xyz))

    def save_displacement_field_back_image(self, filename, split_xyz):
        """Save backward displacement field image to file."""
        assert isinstance(filename, str)
        assert isinstance(split_xyz, bool)
        try_calling(pysirfreg.cSIRFReg_SIRFReg_save_deformation_displacement_image(self.handle, filename, 'back_displacement', split_xyz))

    def update(self):
        """Run the registration"""
        try_calling(pysirfreg.cSIRFReg_SIRFReg_update(self.handle))


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

    def set_transformation_matrix(self, filename):
        """Set transformation matrix."""
        assert isinstance(filename, str)
        _set_char_par(self.handle, self.name, 'transformation_matrix', filename)

    def set_displacement_field(self, displacement_field):
        """Set displacement field."""
        assert isinstance(displacement_field, ImageDataDeformation)
        _setParameter(self.handle, self.name, 'displacement_field', displacement_field.handle)

    def set_deformation_field(self, deformation_field):
        """Set deformation field."""
        assert isinstance(deformation_field, ImageDataDeformation)
        _setParameter(self.handle, self.name, 'deformation_field', deformation_field.handle)

    def set_interpolation_type(self, type):
        """Set interpolation type. 0=nearest neighbour, 1=linear, 3=cubic, 4=sinc."""
        assert isinstance(type, int)
        _set_int_par(self.handle, self.name, 'interpolation_type', type)

    def set_interpolation_type_to_nearestneighbour(self):
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

    def save_resampled_image(self, filename):
        """Save resampled image to file."""
        assert isinstance(filename, str)
        try_calling(pysirfreg.cSIRFReg_SIRFRegNiftyResample_save_resampled_image(self.handle, filename))


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

    def save_image_to_file(self, filename):
        """Save image to file."""
        try_calling(pysirfreg.cSIRFReg_SIRFRegImageWeightedMean3D_save_image_to_file(self.handle, filename))

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

    def save_image_to_file(self, filename):
        """Save image to file."""
        try_calling(pysirfreg.cSIRFReg_SIRFRegImageWeightedMean4D_save_image_to_file(self.handle, filename))

    def get_output(self):
        """Get output."""
        image = ImageData()
        image.handle = _getParameterHandle(self.handle, self.name, 'output')
        check_status(image.handle)
        return image
