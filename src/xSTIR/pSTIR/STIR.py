"""Object-Oriented wrap for the cSTIR-to-Python interface pystir.py."""

# SyneRBI Synergistic Image Reconstruction Framework (SIRF)
# Copyright 2015 - 2022 Rutherford Appleton Laboratory STFC
# Copyright 2015 - 2022 University College London
# Copyright 2019 University of Hull
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
import inspect
import numpy
try:
    import pylab
    HAVE_PYLAB = True
except:
    HAVE_PYLAB = False
import sys
from numbers import Integral, Number
from deprecation import deprecated

from sirf.Utilities import show_2D_array, show_3D_array, error, check_status, \
     try_calling, assert_validity, assert_validities, \
     cpp_int_dtype, cpp_int_array, \
     examples_data_path, existing_filepath, pTest
from sirf import SIRF
from sirf.SIRF import ContiguousError, DataContainer
import sirf.pyiutilities as pyiutil
import sirf.pystir as pystir

import sirf.STIR_params as parms
from sirf.config import SIRF_HAS_NiftyPET
from sirf.config import SIRF_HAS_Parallelproj
from sirf.config import STIR_WITH_CUDA

if sys.version_info[0] >= 3 and sys.version_info[1] >= 4:
    ABC = abc.ABC
else:
    ABC = abc.ABCMeta('ABC', (), {})

INFO_CHANNEL = 0
WARNING_CHANNEL = 1
ERROR_CHANNEL = 2
ALL_CHANNELS = -1

MAX_ACQ_DIMS = 10
MAX_IMG_DIMS = 10


def get_STIR_version_string():
    """Returns STIR engine version as Python str."""
    handle = pystir.cSTIR_STIR_version_string()
    check_status(handle)
    version = pyiutil.charDataFromHandle(handle)
    pyiutil.deleteDataHandle(handle)
    return version


def get_STIR_doc_dir():
    """Returns STIR engine documentation folder name as Python str."""
    handle = pystir.cSTIR_get_STIR_doc_dir()
    check_status(handle)
    path = pyiutil.charDataFromHandle(handle)
    pyiutil.deleteDataHandle(handle)
    return path


def get_STIR_examples_dir():
    """Returns STIR engine examples folder name as Python str."""
    handle = pystir.cSTIR_get_STIR_examples_dir()
    check_status(handle)
    path = pyiutil.charDataFromHandle(handle)
    pyiutil.deleteDataHandle(handle)
    return path


def get_engine_version_string():
    """Returns engine version as Python str."""
    return get_STIR_version_string()


def get_engine_doc_dir():
    """Returns STIR engine documentation folder name as Python str."""
    return get_STIR_doc_dir()


def get_engine_examples_dir():
    """Returns STIR engine examples folder name as Python str."""
    return get_STIR_examples_dir()


def set_verbosity(verbosity):
    """Set the verbosity of all STIR output."""
    try_calling(pystir.cSTIR_setVerbosity(verbosity))


def get_verbosity():
    """Get the verbosity of all STIR output."""
    h = pystir.cSTIR_getVerbosity()
    check_status(h, inspect.stack()[1])
    value = pyiutil.intDataFromHandle(h)
    pyiutil.deleteDataHandle(h)
    return value


def scanner_names():
    h = pystir.cSTIR_scannerNames()
    check_status(h, inspect.stack()[1])
    value = pyiutil.charDataFromHandle(h)
    pyiutil.deleteDataHandle(h)
    return value


def set_max_omp_threads(threads):
    """Set the maximum number of OpenMP threads."""
    try_calling(pystir.cSTIR_setOMPThreads(int(threads)))


def get_max_omp_threads():
    """Get the maximum number of OpenMP threads."""
    h = pystir.cSTIR_getOMPThreads()
    check_status(h, inspect.stack()[1])
    value = pyiutil.intDataFromHandle(h)
    pyiutil.deleteDataHandle(h)
    return value


def set_default_num_omp_threads():
    """Use the default number of OpenMP threads."""
    try_calling(pystir.cSTIR_useDefaultOMPThreads())


def get_default_num_omp_threads():
    """Get default num of OpenMP threads."""
    h = pystir.cSTIR_getDefaultOMPThreads()
    check_status(h, inspect.stack()[1])
    value = pyiutil.intDataFromHandle(h)
    pyiutil.deleteDataHandle(h)
    return value


class MessageRedirector(object):
    """Class for STIR printing redirection to files/stdout/stderr."""

    def __init__(self, info=None, warn='stdout', errr='stderr'):
        """Create MessageRedirector object that redirects STIR's ouput.

        Output produced by info(), warning() and error(0 functions to
        destinations specified respectively by info, warn and err arguments.
        The argument values other than None, stdout, stderr, cout and cerr
        are interpreted as filenames.
        None and empty string value suppresses printing.
        """
        if info is None:
            info = ''
        if not isinstance(info, str):
            raise error(
                'wrong info argument for MessageRedirector constructor')
        elif info in {'stdout', 'stderr', 'cout', 'cerr'}:
            self.info = pystir.newTextPrinter(info)
            self.info_case = 0
        else:
            self.info = pystir.newTextWriter(info)
            self.info_case = 1
        pystir.openChannel(0, self.info)

        if warn is None:
            warn = ''
        if not isinstance(warn, str):
            raise error(
                'wrong warn argument for MessageRedirector constructor')
        elif warn in {'stdout', 'stderr', 'cout', 'cerr'}:
            self.warn = pystir.newTextPrinter(warn)
            self.warn_case = 0
        else:
            self.warn = pystir.newTextWriter(warn)
            self.warn_case = 1
        pystir.openChannel(1, self.warn)

        if errr is None:
            errr = ''
        if not isinstance(errr, str):
            raise error(
                'wrong errr argument for MessageRedirector constructor')
        elif errr in {'stdout', 'stderr', 'cout', 'cerr'}:
            self.errr = pystir.newTextPrinter(errr)
            self.errr_case = 0
        else:
            self.errr = pystir.newTextWriter(errr)
            self.errr_case = 1
        pystir.openChannel(2, self.errr)

    def __del__(self):
        """del."""
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


class Shape(object):
    """Base class for an abstract geometric shape.

    Used as a building block for creating phantom images.
    """

    def __init__(self):
        """init."""
        self.handle = None

    def __del__(self):
        """del."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_origin(self, origin):
        """Sets origin.

        Sets the (discrete) coordinates of
        the shape centre on a voxel grid.
        """
        parms.set_float_par(self.handle, 'Shape', 'x', origin[2])
        parms.set_float_par(self.handle, 'Shape', 'y', origin[1])
        parms.set_float_par(self.handle, 'Shape', 'z', origin[0])

    def get_origin(self):
        """Returns the coordinates of the shape centre on a voxel grid."""
        x = parms.float_par(self.handle, 'Shape', 'x')
        y = parms.float_par(self.handle, 'Shape', 'y')
        z = parms.float_par(self.handle, 'Shape', 'z')
        return (x, y, z)

class Box3D(Shape):
    """Class for elliptic cylinder shape."""

    def __init__(self):
        """init."""
        self.handle = None
        self.name = 'Box3D'
        self.handle = pystir.cSTIR_newObject(self.name)
        check_status(self.handle)

    def __del__(self):
        """del."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_length_x(self, value):
        """Sets dimension x length in mm."""
        parms.set_float_par(self.handle, self.name, 'length_x', value)

    def get_length_x(self):
        """Returns dimension x length in mm."""
        return parms.float_par(self.handle, self.name, 'length_x')

    def set_length_y(self, value):
        """Sets dimension y length in mm."""
        parms.set_float_par(self.handle, self.name, 'length_y', value)

    def get_length_y(self):
        """Returns dimension y length in mm."""
        return parms.float_par(self.handle, self.name, 'length_y')

    def set_length_z(self, value):
        """Sets dimension z length in mm."""
        parms.set_float_par(self.handle, self.name, 'length_z', value)

    def get_length_z(self):
        """Returns dimension z length in mm."""
        return parms.float_par(self.handle, self.name, 'length_z')

    def set_lengths(self, value):
        """Sets lengths in mm."""
        parms.set_float_par(self.handle, self.name, 'length_x', value[0])
        parms.set_float_par(self.handle, self.name, 'length_y', value[1])
        parms.set_float_par(self.handle, self.name, 'length_z', value[2])

    def get_lengths(self):
        """Returns lengths in mm."""
        length_x = parms.float_par(self.handle, self.name, 'length_x')
        length_y = parms.float_par(self.handle, self.name, 'length_y')
        length_z = parms.float_par(self.handle, self.name, 'length_z')
        return (length_x, length_y, length_z)


class Ellipsoid(Shape):
    """Class for ellipsoid shape."""

    def __init__(self):
        """init."""
        self.handle = None
        self.name = 'Ellipsoid'
        self.handle = pystir.cSTIR_newObject(self.name)
        check_status(self.handle)

    def __del__(self):
        """del."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_radius_x(self, value):
        """Sets x radius in mm."""
        parms.set_float_par(self.handle, self.name, 'radius_x', value)

    def get_radius_x(self):
        """Returns x radius in mm."""
        return parms.float_par(self.handle, self.name, 'radius_x')

    def set_radius_y(self, value):
        """Sets y radius in mm."""
        parms.set_float_par(self.handle, self.name, 'radius_y', value)

    def get_radius_y(self):
        """Returns y radius in mm."""
        return parms.float_par(self.handle, self.name, 'radius_y')

    def set_radius_z(self, value):
        """Sets z radius in mm."""
        parms.set_float_par(self.handle, self.name, 'radius_z', value)

    def get_radius_z(self):
        """Returns z radius in mm."""
        return parms.float_par(self.handle, self.name, 'radius_z')


class EllipticCylinder(Shape):
    """Class for elliptic cylinder shape."""

    def __init__(self):
        """init."""
        self.handle = None
        self.name = 'EllipsoidalCylinder'
        self.handle = pystir.cSTIR_newObject(self.name)
        check_status(self.handle)

    def __del__(self):
        """del."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_length(self, value):
        """Sets length in mm."""
        parms.set_float_par(self.handle, self.name, 'length', value)

    def get_length(self):
        """Returns length in mm."""
        return parms.float_par(self.handle, self.name, 'length')

    def set_radius_x(self, value):
        """Sets x radius in mm."""
        parms.set_float_par(self.handle, self.name, 'radius_x', value)

    def get_radius_x(self):
        """Returns x radius in mm."""
        return parms.float_par(self.handle, self.name, 'radius_x')

    def set_radius_y(self, value):
        """Sets y radius in mm."""
        parms.set_float_par(self.handle, self.name, 'radius_y', value)

    def get_radius_y(self):
        """Returns y radius in mm."""
        return parms.float_par(self.handle, self.name, 'radius_y')

    def set_radii(self, radii):
        """Sets radii in mm."""
        parms.set_float_par(self.handle, self.name, 'radius_x', radii[1])
        parms.set_float_par(self.handle, self.name, 'radius_y', radii[0])

    def get_radii(self):
        """Returns radii in mm."""
        rx = parms.float_par(self.handle, self.name, 'radius_x')
        ry = parms.float_par(self.handle, self.name, 'radius_y')
        return (rx, ry)


# class ImageData(DataContainer):
class ImageData(SIRF.ImageData):
    """Class for PET image data objects.

    ImageData objects contain both geometric data and the actual voxel
    values. You have to use the `as_array` method to get an array with
    the voxel values, and use the `fill` function to change the voxel values.
    """

    def __init__(self, arg=None):
        """Creates an ImageData object.

        Arguments:
            arg : Python str or AcquisitionData or None, interpreted as follows:
        str            : read the object from a file specified by <arg>
                         (the file format has to be support by STIR).
        STIR.AcquisitionData: create an object compatible with the scanner data
                         recorded in an AcquisitionData object <arg>.
                         This sets default voxel sizes.
        SIRF.ImageData : (attempt to) convert data from another SIRF ImageData container.
        None           : create an empty ImageData object. Call initialise()
                         method before using it.
        """
        self.handle = None
        if isinstance(arg, str):
            self.handle = pystir.cSTIR_objectFromFile('Image', arg)
            check_status(self.handle)
        elif isinstance(arg, AcquisitionData):
            if arg.handle is None:
                raise AssertionError()
            self.handle = pystir.cSTIR_imageFromAcquisitionData(arg.handle)
            check_status(self.handle)
        elif isinstance(arg, SIRF.ImageData):
            if arg.handle is None:
                raise AssertionError()
            self.handle = pystir.cSTIR_imageFromImageData(arg.handle)
            check_status(self.handle)
        elif arg is not None:
            raise error(
                'wrong argument ' + repr(arg) + ' for ImageData constructor')
        self.name = 'ImageData'
        self.rimsize = -1

    def __del__(self):
        """Deallocates this ImageData object."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def same_object(self):
        """See DataContainer method."""
        return ImageData()

    def get_info(self):
        """Returns the STIR metadata as Python str."""
        handle = pystir.cSTIR_get_info(self.handle)
        check_status(handle)
        info = pyiutil.charDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return info

    @property
    def modality(self):
        """Returns imaging modality as Python string."""
        return parms.char_par(self.handle, 'ImageData', 'modality')

    def set_modality(self, mod):
        """Sets imaging modality.

        mod: "PT" or "NM" or "MR" or "CT" or "US" or "Optical"
        """
        return parms.set_char_par(self.handle, 'ImageData', 'modality', mod)

    @property
    def __array_interface__(self):
        """As per https://numpy.org/doc/stable/reference/arrays.interface.html"""
        if not self.supports_array_view:
            raise ContiguousError("views not supported, please consider using `asarray()` or `as_array()`")
        return {'shape': self.shape, 'typestr': '<f4', 'version': 3,
                'data': (parms.size_t_par(self.handle, 'ImageData', 'address'), False)}

    def initialise(self, dim, vsize=(1., 1., 1.), origin=(0., 0., 0.)):
        """
        Sets image size and geometric information.

        If self is not an empty ImageData, its contents are erased.
        Dimemsions (number of voxels) are required,
        spacing and offset are optional.
        Spacing and offset are to be given in mm.
        All should be given as (z,y,x).
        """
        if not isinstance(dim, tuple):
            raise error(
                "sirf.STIR.Imagedata.initialise: dim should be a tuple")
        if not isinstance(vsize, tuple):
            raise error(
                "sirf.STIR.Imagedata.initialise: vsize should be a tuple")
        if not isinstance(origin, tuple):
            raise error(
                "sirf.STIR.Imagedata.initialise: origin should be a tuple")
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
        self.handle = None
        voxels = pystir.cSTIR_voxels3DF(
            int(dim[2]), int(dim[1]), int(dim[0]),
            float(vsize[2]), float(vsize[1]), float(vsize[0]),
            float(origin[2]), float(origin[1]), float(origin[0]))
        check_status(voxels)
        self.handle = pystir.cSTIR_imageFromVoxels(voxels)
        check_status(self.handle)
        pyiutil.deleteDataHandle(voxels)

    def fill(self, value):
        """Sets the voxel values.

        The argument is either ImageData or 3D Numpy ndarray of values or a
        scalar to be assigned at each voxel. When using an ndarray, the array
        must have the same size as an array returned by `as_array`.
        """
        if self.handle is None:
            raise AssertionError()
        if isinstance(value, ImageData):
            super(ImageData, self).fill(value)
#            try_calling(pystir.cSTIR_setImageDataFromImage \
#                        (self.handle, value.handle))
        elif isinstance(value, numpy.ndarray):
            dims = self.dimensions()
            shape = value.shape
            if shape != dims:
                msg = 'cannot fill ImageData of size %s' \
                      + ' with data of size %s'
                raise ValueError(msg % (repr(dims), repr(shape)))
            if value.dtype is numpy.dtype('float32'):
                # print('keeping dtype float32')
                v = value
            else:
                # print('changing dtype to float32')
                v = value.astype(numpy.float32)
            if not v.flags['C_CONTIGUOUS']:
                v = numpy.ascontiguousarray(v)
            try_calling(pystir.cSTIR_setImageData(self.handle, v.ctypes.data))
        elif isinstance(value, (Number, numpy.number)):
            try_calling(pystir.cSTIR_fillImage(self.handle, float(value)))
        else:
            raise TypeError('wrong fill value.' + \
                        ' Should be ImageData, numpy.ndarray, float or int. Got {}'\
                        .format(type(value)))
        return self

    def get_uniform_copy(self, value=1.0):
        """Creates a copy of this image filled with <value>."""
        if self.handle is None:
            raise AssertionError()
        image = ImageData()
        image.handle = pystir.cSTIR_imageFromImage(self.handle)
        check_status(image.handle)
        image.fill(value)
        return image

    def add_shape(self, shape, scale, num_samples_in_each_direction=1):
        """Adds a shape to self - see Shape above.

        If a shape partially fills a voxel, it is possible to choose the
        number of samples that will be used in each direction to determine the
        fraction of the voxel that is filled by the shape. For a 3D image,
        using num_samples_in_each_direction=2 would result in 2^3=8 samples.
        """
        if self.handle is None:
            raise AssertionError()
        assert_validity(shape, Shape)
        try_calling(pystir.cSTIR_addShape(
            self.handle, shape.handle, scale,
            int(num_samples_in_each_direction)))

    def read_from_file(self, filename):
        """
        Reads data from file.

        Replaces the current content of the object.
        """
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
        self.handle = pystir.cSTIR_objectFromFile('Image', filename)
        check_status(self.handle)

    def dimensions(self):
        """Returns image dimensions as a tuple (nz, ny, nx)."""
        if self.handle is None:
            raise AssertionError()
        dim = numpy.ndarray((MAX_IMG_DIMS,), dtype=cpp_int_dtype())
        try_calling(
            pystir.cSTIR_getImageDimensions(self.handle, dim.ctypes.data))
        return tuple(dim[:3])  # [::-1])

    def voxel_sizes(self):
        """Returns image voxel sizes as a tuple (vz, vy, vx)."""
        if self.handle is None:
            raise AssertionError()
        vs = numpy.ndarray((3,), dtype=numpy.float32)
        try_calling(
            pystir.cSTIR_getImageVoxelSizes(self.handle, vs.ctypes.data))
        return tuple(vs)  # [::-1])

    @property
    def spacing(self):
        return self.voxel_sizes()

    def transf_matrix(self):
        """Returns transformation matrix."""
        if self.handle is None:
            raise AssertionError()
        tm = numpy.ndarray((4, 4), dtype=numpy.float32)
        try_calling(
            pystir.cSTIR_getImageTransformMatrix(self.handle, tm.ctypes.data))
        return tm

    def as_array(self):
        """
        WARNING: you probably should use `.asarray()` (no underscore) instead.

        Returns 3D Numpy ndarray with values at the voxels.
        """
        if self.handle is None:
            raise AssertionError()
        array = numpy.ndarray(self.dimensions(), dtype=numpy.float32)
        try_calling(pystir.cSTIR_getImageData(self.handle, array.ctypes.data))
        return array

    def write_par(self, filename, par):
        """Writes with parameter file."""
        try_calling(pystir.cSTIR_writeImage_par(self.handle, filename, par))

    def show(self, slice=None, title=None):
        """Displays xy-cross-section(s) of this image."""
        if self.handle is None:
            raise AssertionError()
        if not HAVE_PYLAB:
            print('pylab not found')
            return
        data = self.as_array()
        nz = data.shape[0]
        if isinstance(slice, (Integral,numpy.integer)):
            if slice < 0 or slice >= nz:
                raise IndexError('Slice index out of range')
            show_2D_array('slice %d' % slice, data[slice, :, :])
            return
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
            show_3D_array(data, index=slice[f : t], tile_shape=tiles,
                          label='slice', xlabel='x', ylabel='y',
                          suptitle=title, show=(t == ni))
            f = t

    def zoom_image(self, zooms=(1., 1., 1.), offsets_in_mm=(0., 0., 0.),
                   size=(-1, -1, -1), scaling='preserve_sum'):
        """
        Returns a zoomed image.

        All coordinates and indices are given as (z,y,x).
        To leave the size unchanged in any dimension, set the
        corresponding size to -1
        Supported scaling options are: 'preserve_sum', 'preserve_values' and
        'preserve_projections'
        """
        zoomed_im = self.clone()

        if not isinstance(zooms, tuple):
            raise error('zoom_image: zooms should be tuple')
        if not isinstance(offsets_in_mm, tuple):
            raise error('zoom_image: offsets_in_mm should be tuple')
        if not isinstance(size, tuple):
            raise error('zoom_image: size should be tuple')
        np_zooms = numpy.asarray(zooms, dtype=numpy.float32)
        np_offsets_in_mm = numpy.asarray(offsets_in_mm, dtype=numpy.float32)
        np_size = numpy.asarray(size, dtype=cpp_int_dtype())

        try_calling(pystir.cSTIR_ImageData_zoom_image(
            zoomed_im.handle, np_zooms.ctypes.data,
            np_offsets_in_mm.ctypes.data, np_size.ctypes.data, scaling))

        return zoomed_im

    def zoom_image_as_template(self, template_image, scaling='preserve_sum'):
        """
        Returns a zoomed image based on a template image's geometry.

        Supported scaling options are: 'preserve_sum', 'preserve_values' and
        'preserve_projections'
        """
        zoomed_image = template_image.clone()

        if not isinstance(template_image, ImageData):
            raise error('zoom_image_as_template: template should be ImageData')

        ### because of a bug somewherre
        try_calling(pystir.cSTIR_ImageData_zoom_image_as_template(
            zoomed_image.handle, self.handle, scaling))

        return zoomed_image

    def move_to_scanner_centre(self, proj_data):
        """Moves the image to the scanner centre.

        AcquisitionData is required as bed shift etc. will be taken into
        account when available.
        """
        if not isinstance(proj_data, AcquisitionData):
            raise error('move_to_scanner_centre: proj_data should' +
                        'be of type AcquisitionData')
        if proj_data.handle is None:
            raise error('move_to_scanner_centre: proj_data is not initialised')
        if self.handle is None:
            raise error('move_to_scanner_centre: image is not initialised')
        moved_im = self.clone()
        try_calling(pystir.cSTIR_ImageData_move_to_scanner_centre(
            moved_im.handle, proj_data.handle))

        return moved_im

    @property
    def shape(self):
        return self.dimensions()


SIRF.ImageData.register(ImageData)


class ImageDataProcessor(object):
    """Class for image processors.

    An ImageDataProcessor changes an image in some way, e.g. by filtering.
    """

    def __init__(self):
        """init."""
        self.handle = None
        # TODO: handle input and output in cSTIR
        self.input = None
        self.output = None

    def __del__(self):
        """del."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def apply(self, image):
        """Applies this filter to the specified image."""
        assert_validity(image, ImageData)
        try_calling(pystir.cSTIR_applyImageDataProcessor(
            self.handle, image.handle))

    def set_input(self, input):
        """Sets the input data."""
        assert_validity(input, ImageData)
        self.input = input

    def process(self, input=None):
        """Processes data."""
        if input is not None:
            self.input = input
        if self.input is None:
            raise error('input image not set')
        assert_validity(self.input, ImageData)
        self.output = self.input.clone()
        self.apply(self.output)
        return self.output

    def get_output(self):
        """Returns the output data."""
        return self.output


class SeparableGaussianImageFilter(ImageDataProcessor):
    """Implements Gaussian filtering.

    The filtering operation is performed as 3 separate one-dimensional filters
    in each spacial direction.
    """

    def __init__(self):
        """init."""
        self.handle = None
        self.input = None
        self.output = None
        self.name = 'SeparableGaussianImageFilter'
        self.handle = pystir.cSTIR_newObject(self.name)
        check_status(self.handle)

    def __del__(self):
        """del."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_up(self, image):
        """Sets up."""
        assert_validity(image, ImageData)
        try_calling(pystir.cSTIR_setupImageDataProcessor(
            self.handle, image.handle))

    def set_fwhms(self, fwhms):
        """Sets FWHM in mm."""
        parms.set_float_par(self.handle, self.name, 'fwhms_x', fwhms[2])
        parms.set_float_par(self.handle, self.name, 'fwhms_y', fwhms[1])
        parms.set_float_par(self.handle, self.name, 'fwhms_z', fwhms[0])

    def set_max_kernel_sizes(self, mks):
        """Sets max kernel sizes in voxels."""
        parms.set_int_par(self.handle, self.name, 'max_kernel_size_x', mks[2])
        parms.set_int_par(self.handle, self.name, 'max_kernel_size_y', mks[1])
        parms.set_int_par(self.handle, self.name, 'max_kernel_size_z', mks[0])

    def set_normalise(self, norm=True):
        """Sets normalise."""
        v = 1 if norm else 0
        parms.set_int_par(self.handle, self.name, 'normalise', v)


class TruncateToCylinderProcessor(ImageDataProcessor):
    """Class for the image filter that zeroes the image outside the cylinder
       of the same xy-diameter and z-size as those of the image.
    """

    def __init__(self):
        """init."""
        self.handle = None
        self.name = 'TruncateToCylindricalFOVImageProcessor'
        self.handle = pystir.cSTIR_newObject(self.name)
        check_status(self.handle)

    def __del__(self):
        """del."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_strictly_less_than_radius(self, flag):
        """Specifies the behaviour on the boundary.

        Specifies whether the area not affected by filtering is strictly
        inside the cylinder (flag = True) or not (flag = False).
        """
        parms.set_char_par(
            self.handle, 'TruncateToCylindricalFOVImageProcessor',
            'strictly_less_than_radius', repr(flag))

    def get_strictly_less_than_radius(self):
        """Returns the boundary behaviour flag.

        Returns the answer to the question: Is the area not affected by
        filtering strictly inside the cylinder?
        """
        return parms.int_par(
            self.handle, 'TruncateToCylindricalFOVImageProcessor',
            'strictly_less_than_radius') != 0


class RayTracingMatrix(object):
    """Sparse matrix ray-tracting projector.

    Class for objects holding sparse matrix representation of the ray
    tracing projector G (see AcquisitionModel class).
    """

    name = 'RayTracingMatrix'

    def __init__(self):
        """init."""
        self.handle = pystir.cSTIR_newObject(self.name)
        check_status(self.handle)
        parms.set_int_par(self.handle, self.name, 'num_tangential_LORs', 2)

    def __del__(self):
        """del."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def get_info(self):
        """Returns the metadata from STIR as Python str."""
        handle = pystir.cSTIR_get_MatrixInfo(self.handle)
        check_status(handle)
        info = pyiutil.charDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return info

    def set_num_tangential_LORs(self, value):
        """Sets the number of tangential LORs.

        Sets the number of LORs (or rays) for each bin in the sinogram.
        They are currently (approximately) parallel and spaced in the
        tangential direction (i.e. orthogonal to the axial direction).
        """
        parms.set_int_par(self.handle, self.name, 'num_tangential_LORs', value)
        return self

    def get_num_tangential_LORs(self):
        """Returns the number of LORs for each bin in the sinogram."""
        return parms.int_par(self.handle, self.name, 'num_tangential_LORs')

    def enable_cache(self, value=True):
        """Enables or disables the caching mechanism."""
        parms.set_bool_par(self.handle, self.name, 'enable_cache', value)
        return self

    def set_restrict_to_cylindrical_FOV(self, value=True):
        """Enables or disables using a circular axial FOV (vs rectangular)."""
        parms.set_bool_par(self.handle, self.name, 'restrict_to_cylindrical_FOV', value)
        return self

    def set_do_symmetry_90degrees_min_phi(self, value=True):
        """Enables or disables a symmetry (disabling saves memory but might increase computation time)."""
        parms.set_bool_par(self.handle, self.name, 'do_symmetry_90degrees_min_phi', value)
        return self

    def set_do_symmetry_180degrees_min_phi(self, value=True):
        """Enables or disables a symmetry (disabling saves memory but might increase computation time)."""
        parms.set_bool_par(self.handle, self.name, 'do_symmetry_180degrees_min_phi', value)
        return self

    def set_do_symmetry_swap_segment(self, value=True):
        """Enables or disables a symmetry (disabling saves memory but might increase computation time)."""
        parms.set_bool_par(self.handle, self.name, 'do_symmetry_swap_segment', value)
        return self

    def set_do_symmetry_swap_s(self, value=True):
        """Enables or disables a symmetry (disabling saves memory but might increase computation time)."""
        parms.set_bool_par(self.handle, self.name, 'do_symmetry_swap_s', value)
        return self

    def set_do_symmetry_shift_z(self, value=True):
        """Enables or disables a symmetry (disabling saves memory but might increase computation time)."""
        parms.set_bool_par(self.handle, self.name, 'do_symmetry_shift_z', value)
        return self

class SPECTUBMatrix:
    '''
    Class for objects holding sparse matrix representation of a SPECT
    projector (developed at the University of Barcelona) (see AcquisitionModel class).
    '''
    name = 'SPECTUBMatrix'

    def __init__(self):
        '''
        Create a new matrix. Default settings use neither attenuation nor resolution modelling.
        '''
        self.handle = pystir.cSTIR_newObject(self.name)
        check_status(self.handle)

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def set_keep_all_views_in_cache(self, value):
        '''
        Enable keeping the matrix in memory.

        This speeds-up the calculations, but can use a lot of memory.

        You have to call set_up() after this (unless the value didn't change).
        '''
        parms.set_int_par(self.handle, self.name, 'keep_all_views_in_cache', value)
        return self
    def get_keep_all_views_in_cache(self):
        '''
        Returns a bool checking if we're keeping the whole matrix in memory or not.
        '''
        return parms.int_par(self.handle, self.name, 'keep_all_views_in_cache') != 0
    def set_attenuation_image(self, value):
        '''
        Sets the attenuation image used by the projector.
        '''
        assert_validity(value, ImageData)
        parms.set_parameter(self.handle, self.name, 'attenuation_image', value.handle)
        return self
    def get_attenuation_image(self):
        '''
        Returns the attenuation image used by the projector.
        '''
        image = ImageData()
        image.handle = parms.parameter_handle(self.handle, self.name, 'attenuation_image')
        return image

    def set_resolution_model(self, collimator_sigma_0_in_mm, collimator_slope_in_mm, full_3D = True):
        '''
        Set the parameters for the depth-dependent resolution model

        The detector and collimator blurring is modelled as a Gaussian with sigma dependent on the
        distance from the collimator.

        sigma_at_depth = collimator_slope * depth_in_mm + collimator sigma 0

        Set slope and sigma_0 to zero to avoid resolution modelling.

        You have to call set_up() after this.
        '''
        try_calling(pystir.cSTIR_SPECTUBMatrixSetResolution(self.handle, collimator_sigma_0_in_mm, collimator_slope_in_mm, full_3D))

    def set_up(self, acq, img):
        try_calling(pystir.cSTIR_setupSPECTUBMatrix(self.handle, acq.handle, img.handle))


class PinholeSPECTUBMatrix:
    """
    Class for objects holding sparse matrix representation of a pinhole SPECT
    projector (developed at the University of Barcelona) (see AcquisitionModel class).
    """
    name = 'PinholeSPECTUBMatrix'

    def __init__(self):
        """Create a new matrix. Default settings use neither attenuation, PSF, or DOI modelling."""
        self.handle = pystir.cSTIR_newObject(self.name)
        check_status(self.handle)

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def get_maximum_number_of_sigmas(self):
        """Returns the number of sigmas to consider when correcting for intrinsic PSF."""
        return parms.float_par(self.handle, self.name, 'maximum_number_of_sigmas')

    def set_maximum_number_of_sigmas(self, value):
        """Sets the number of sigmas to consider when correcting for intrinsic PSF."""
        parms.set_float_par(self.handle, self.name, 'maximum_number_of_sigmas', value)

    def get_spatial_resolution_PSF(self):
        """Returns the spatial high resolution in which to sample distributions (in cm)."""
        return parms.float_par(self.handle, self.name, 'spatial_resolution_PSF')

    def set_spatial_resolution_PSF(self, value):
        """Sets the spatial high resolution in which to sample distributions (in cm)."""
        parms.set_float_par(self.handle, self.name, 'spatial_resolution_PSF', value)

    def get_subsampling_factor_PSF(self):
        """Returns the subsampling factor to compute convolutions when PSF or DOI corrections are enabled."""
        return parms.int_par(self.handle, self.name, 'subsampling_factor_PSF')

    def set_subsampling_factor_PSF(self, value):
        """Sets the subsampling factor to compute convolutions when PSF or DOI corrections are enabled."""
        parms.set_int_par(self.handle, self.name, 'subsampling_factor_PSF', value)

    def set_detector_file(self, filename):
        """Sets the name of the file containing the detector information."""
        parms.set_char_par(self.handle, self.name, 'detector_file', filename)

    def set_collimator_file(self, filename):
        """Sets the name of the file containing the collimator information."""
        parms.set_char_par(self.handle, self.name, 'collimator_file', filename)

    def get_psf_correction(self):
        """Returns the setting for enabling corrections for intrinsic PSF."""
        return parms.char_par(self.handle, self.name, 'psf_correction')

    def set_psf_correction(self, value):
        """Enable or disable corrections for intrinsic PSF."""
        parms.set_char_par(self.handle, self.name, 'psf_correction', value)

    def get_doi_correction(self):
        """Returns the setting for enabling corrections for depth of interaction."""
        return parms.char_par(self.handle, self.name, 'doi_correction')

    def set_doi_correction(self, value):
        """Enable or disable corrections for depth of interaction."""
        parms.set_char_par(self.handle, self.name, 'doi_correction', value)

    def get_attenuation_type(self):
        """Returns the attenuation type: full, simple, or no."""
        return parms.char_par(self.handle, self.name, 'attenuation_type')

    def set_attenuation_type(self, value):
        """Set the attenuation type to full, simple, or no."""
        parms.set_char_par(self.handle, self.name, 'attenuation_type', value)

    def get_attenuation_image(self):
        """Returns the attenuation image used by the projector."""
        image = ImageData()
        image.handle = parms.parameter_handle(self.handle, self.name, 'attenuation_image')
        return image

    def set_attenuation_image(self, value):
        """Sets the attenuation image used by the projector."""
        assert_validity(value, ImageData)
        parms.set_parameter(self.handle, self.name, 'attenuation_image', value.handle)
        return self

    def get_object_radius(self):
        """Returns the radius of the object in the xy plane of the image volume."""
        return parms.float_par(self.handle, self.name, 'object_radius')

    def set_object_radius(self, value):
        """Sets the radius of the object in the xy plane of the image volume. Could be used for masking."""
        parms.set_float_par(self.handle, self.name, 'object_radius', value)

    def get_mask_image(self):
        """Returns the mask image used by the projector."""
        image = ImageData()
        image.handle = parms.parameter_handle(self.handle, self.name, 'mask_image')
        return image

    def set_mask_image(self, value):
        """Sets the mask image used by the projector."""
        assert_validity(value, ImageData)
        parms.set_parameter(self.handle, self.name, 'mask_image', value.handle)
        return self

    def get_keep_all_views_in_cache(self):
        """Returns a bool checking if we're keeping the whole matrix in memory or not."""
        return parms.bool_par(self.handle, self.name, 'keep_all_views_in_cache')

    def set_keep_all_views_in_cache(self, value):
        """Enable keeping the matrix in memory to speed-up calculations (can use lots of memory)."""
        parms.set_bool_par(self.handle, self.name, 'keep_all_views_in_cache', value)
        return self

    def get_mask_from_attenuation_map(self):
        """Returns a bool checking if we're masking with the attenuation map or not."""
        return parms.bool_par(self.handle, self.name, 'mask_from_attenuation_map')

    def set_mask_from_attenuation_map(self, value):
        """Enable masking from attenuation map if mask file is not set."""
        parms.set_bool_par(self.handle, self.name, 'mask_from_attenuation_map', value)
        return self


class ScanData(DataContainer):
    """Abstract base class for STIR raw data."""

    def get_info(self):
        """Returns the STIR metadata as Python str."""
        handle = pystir.cSTIR_get_info(self.handle)
        check_status(handle)
        info = pyiutil.charDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return info


class ListmodeData(ScanData):
    """Class for STIR list mode data."""

    def __init__(self, filename=None):
        self.handle = None
        self.name = 'ListmodeData'
        self.read_only = False
        if filename is None:
            return
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
        self.handle = pystir.cSTIR_objectFromFile('ListmodeData', filename)
        check_status(self.handle)
        self.read_only = True

    def __del__(self):
        """del."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def read_from_file(self, filename):
        """
        Reads data from file.

        Replaces the current content of the object.
        """
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
        self.handle = pystir.cSTIR_objectFromFile('ListmodeData', filename)
        check_status(self.handle)
        self.read_only = True

    def acquisition_data_template(self):
        """
        Construct an AcquisitionData object corresponding to the listmode data
        if no additional compression (such as as mashing, or rebinning) is used.
        """
        if self.handle is None:
            raise AssertionError('ListmodeData not yet set')
        acq_data = AcquisitionData()
        acq_data.handle = pystir.cSTIR_acquisitionDataFromListmode(self.handle)
        return acq_data


ScanData.register(ListmodeData)


class AcquisitionData(ScanData):
    """Class for stir sinogram (or projection) acquisition data."""

    def __init__(self, src=None, span=1, max_ring_diff=-1, view_mash_factor=1, tof_mash_factor=1):
        """Creates new AcquisitionData.

        Can create object from a file or another AcquisitionData object.
        src:  file name (Python str) or AcquisitionData object or scanner name
        """
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
                self.handle = pystir.cSTIR_objectFromFile(
                    'AcquisitionData', src)
                self.read_only = self.get_storage_scheme() == 'file'
                self.src = 'file'
            else:
                # src is a scanner name
                self.handle = pystir.cSTIR_acquisitionDataFromScannerInfo(
                    src, span, max_ring_diff, view_mash_factor, tof_mash_factor)
                if pyiutil.executionStatus(self.handle) != 0:
                    msg = pyiutil.executionError(self.handle)
                    if msg == 'Unknown scanner':
                        raise error(
                            'Unknown scanner ' + src +
                            ' or missing raw data file extension')
                self.src = 'scanner'
        elif isinstance(src, AcquisitionData):
            # src is AcquisitionData
            if src.handle is None:
                raise AssertionError()
            self.handle = pystir.cSTIR_acquisitionDataFromTemplate(src.handle)
            self.src = 'template'
        else:
            raise error('Wrong source in AcquisitionData constructor')
        check_status(self.handle)

    def __del__(self):
        """del."""
        # print('deleting AcquisitionData object originated from ', self.src)
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    @staticmethod
    def set_storage_scheme(scheme):
        """Sets acquisition data storage scheme.

        scheme = 'file' (default):
            all acquisition data generated from now on will be kept in
            scratch files deleted after the user's script terminates
        scheme = 'memory':
            all acquisition data generated from now on will be kept in RAM
            (avoid if data is very large)
        """
        try_calling(pystir.cSTIR_setAcquisitionDataStorageScheme(scheme))

    @staticmethod
    def get_storage_scheme():
        """Returns acquisition data storage scheme."""
        handle = pystir.cSTIR_getAcquisitionDataStorageScheme()
        check_status(handle)
        scheme = pyiutil.charDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return scheme

    def same_object(self):
        """See DataContainer method."""
        return AcquisitionData()

    def read_from_file(self, filename):  # 'read_from_file' is misleading
        """
        Reads data from file.

        Replaces the current content of the object.
        """
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
        self.handle = pystir.cSTIR_objectFromFile('AcquisitionData', filename)
        check_status(self.handle)
        self.read_only = True

    def create_uniform_image(self, value=0, xy=None):
        """Crates uniform image.

        Creates ImageData object containing PET image of z-dimension
        and voxel sizes compatible with the scanner geometry stored
        in this AcquisitionData object and assigns a given value
        to all voxels;
        value: a Python float.
        xy   : x and y dimensions tuple (if None, set by STIR)
        """
        if self.handle is None:
            raise AssertionError()
        image = ImageData()
        if xy is None:
            image.handle = pystir.cSTIR_imageFromAcquisitionData(self.handle)
        elif isinstance(xy, tuple):
            image.handle = pystir.cSTIR_imageFromAcquisitionDataAndNxNy(
                self.handle, xy[1], xy[0])
        elif isinstance(xy, int):
            image.handle = pystir.cSTIR_imageFromAcquisitionDataAndNxNy(
                self.handle, xy, xy)
        else:
            raise error('Wrong second argument in create_uniform_image')
        check_status(image.handle)
        image.fill(value)
        return image

    def dimensions(self):
        """Returns a tuple of the data dimensions.

        Contains:
        - number of TOF bins
        - number of sinograms
        - number of views
        - number of tangential positions.
        """
        if self.handle is None:
            raise AssertionError()
        dim = numpy.ndarray((MAX_ACQ_DIMS,), dtype=cpp_int_dtype())
        try_calling(pystir.cSTIR_getAcquisitionDataDimensions(
            self.handle, dim.ctypes.data))
        dim = dim[:4]
        return tuple(dim[::-1])

    def get_tof_mash_factor(self):
        '''Returns TOF mashing factor.'''
        return parms.int_par(self.handle, 'AcquisitionData', 'tof_mash_factor')

    @property
    def modality(self):
        """Returns imaging modality as Python string."""
        return parms.char_par(self.handle, 'AcquisitionData', 'modality')

    @property
    def __array_interface__(self):
        """As per https://numpy.org/doc/stable/reference/arrays.interface.html"""
        if not self.supports_array_view:
            raise ContiguousError("views not supported, please consider using `asarray()` or `as_array()`")
        return {'shape': self.shape, 'typestr': '<f4', 'version': 3,
                'data': (parms.size_t_par(self.handle, 'AcquisitionData', 'address'), False)}

    def as_array(self):
        """
        WARNING: you probably should use `.asarray()` (no underscore) instead.

        Returns bin values as ndarray.

        Return a copy of acquisition data stored in this object as a
        NumPy ndarray of 4 dimensions (in default C ordering of data):
        - number of TOF bins
        - number of sinograms
        - number of views
        - number of tangential positions.
        """
        if self.handle is None:
            raise AssertionError()
        array = numpy.ndarray(self.dimensions(), dtype=numpy.float32)
        try_calling(pystir.cSTIR_getAcquisitionData(
            self.handle, array.ctypes.data))
        return array

    def fill(self, value):
        """Fills the object with values.

        value:  either NumPy ndarray or another AcquisitionData object
                or Python float.
        """
        if self.handle is None:
            raise AssertionError()
        if self.read_only:
            raise error(
                'Cannot fill read-only object, consider filling a clone')
        if isinstance(value, numpy.ndarray):
            dims = self.dimensions()
            shape = value.shape
            if shape != dims:
                msg = 'cannot fill AcquisitionData of size %s' \
                      + ' with data of size %s'
                raise ValueError(msg % (repr(dims), repr(shape)))
            if value.dtype is numpy.dtype('float32'):
                # print('keeping dtype float32')
                v = value
            else:
                # print('changing dtype to float32')
                v = value.astype(numpy.float32)
            if not v.flags['C_CONTIGUOUS']:
                v = numpy.ascontiguousarray(v)
            try_calling(pystir.cSTIR_setAcquisitionData(
                self.handle, v.ctypes.data))
        elif isinstance(value, AcquisitionData):
            if value.handle is None:
                raise AssertionError()
            try_calling(pystir.cSTIR_fillAcquisitionDataFromAcquisitionData(
                self.handle, value.handle))
        elif isinstance(value, float):
            try_calling(pystir.cSTIR_fillAcquisitionData(self.handle, value))
        elif isinstance(value, (Integral,numpy.number)):
            try_calling(pystir.cSTIR_fillAcquisitionData(
                self.handle, float(value)))
        else:
            raise TypeError('Wrong fill value.' + \
                ' Should be numpy.ndarray, AcquisitionData, float or int, got {}'\
                .format(type(value)))
        return self

    def get_uniform_copy(self, value=0):
        """Returns a copy of this object filled with given value.

        Returns a true copy of this object filled with a given value;
        value:  a Python float.
        """
        ad = AcquisitionData(self)
        ad.fill(value)
        ad.src = 'copy'
        return ad

    def rebin(self, num_segments_to_combine,
              num_views_to_combine=1, num_tang_poss_to_trim=0,
              do_normalisation=True, max_in_segment_num_to_process=-1,
              num_tof_bins_to_combine=1):
        """Re-bins the data to lower resolution.

        Keyword arguments:
		num_segments_to_combine -- combines multiple oblique 'segments' together. If set to the
		    total number of segments, this corresponds to SSRB. Another example is if the input data
			has 'span=1', the output span will be equal to the num_segments_to_combine.
		num_views_to_combine -- combines neighbouring views. Needs to be a divisor of the total
		    number of views in the data.
		num_tang_poss_to_trim -- removes a number of tangential positions (horizontal direction
		    in the sinogram) at each end
		do_normalisation -- if True, averages the data, otherwise it adds the data. Often
		    the latter is required for emission data (as it preserves Poisson statistics),
			while the former should be used for corrected data (or for attenuation correction factors).
		max_in_segment_num_to_process -- by default all input data are used. If set to a non-negative
		    number, it will remove the most oblique segments.
		num_tof_bins_to_combine -- number of TOF bins to combine.
        """
        ad = AcquisitionData()
        ad.handle = pystir.cSTIR_rebinnedAcquisitionData(
            self.handle,
            num_segments_to_combine, num_views_to_combine,
            num_tang_poss_to_trim, do_normalisation,
            max_in_segment_num_to_process, num_tof_bins_to_combine)
        check_status(ad.handle)
        return ad

    def show(self, sino=None, tof=0, title=None):
        '''Displays selected sinograms.'''
        if self.handle is None:
            raise AssertionError()
        if not HAVE_PYLAB:
            print('pylab not found')
            return
        data = self.as_array()
        if tof <0 or tof >= data.shape[0]:
            raise IndexError('TOF bin index out of range')
        nz = data.shape[1]
        if isinstance(sino, (Integral,numpy.integer)):
            if sino < 0 or sino >= nz:
                raise IndexError('Slice index out of range')
            show_2D_array('sinogram %d' % sino, data[tof, sino, :, :])
            return
        elif sino is None:
            ns = nz
            sino = range(nz)
        else:
            try:
                ns = len(sino)
            except:
                raise error('wrong sinograms list')
        if title is None:
            title = 'Selected sinograms'
        if ns >= 16:
            tiles = (4, 4)
        else:
            tiles = None
        f = 0
        while f < ns:
            t = min(f + 16, ns)
            show_3D_array(
                data[0, :, :, :],
                index=sino[f: t], tile_shape=tiles,
                label='sinogram',
                xlabel='tang.pos', ylabel='view',
                suptitle=title, show=(t == ns))
            f = t

    def get_subset(self, views):
        """Returns the subset of self data formed by specified views

        views: array of views (will be converted to numpy ndarray)
        """
        # Ensure the array passed to C++ is a contiguous array of C++ int's
        v = cpp_int_array(views)
        n = len(views)
        subset = AcquisitionData()
        subset.handle = pystir.cSTIR_get_subset(self.handle, n, v.ctypes.data)
        check_status(subset.handle)
        return subset

    @property
    def shape(self):
        return self.dimensions()


ScanData.register(AcquisitionData)


class ListmodeToSinograms(object):
    """
    Class for listmode-to-sinogram converter.

    This class reads list mode data and produces corresponding *sinograms*,
    i.e. histogrammed data in the format of PETAcquisitionData.

    It has two main functions:
      - process() can be used to read prompts and/or delayed coincidences to
        produce a single PETAcquisitionData.
        Two conversion flags decide what is to be done with 3 possible cases:
        - `store_prompts`=`true`, `store_delayeds`=`false`:
        only prompts stored
        - `store_prompts`=`false`, `store_delayeds`=`true`:
        only delayeds stored
        - `store_prompts`=`true`, `store_delayeds`=`true`:
        prompts-delayeds stored
        Clearly, enabling the `store_delayeds` option only makes sense if the
        data was acquired accordingly.
      - estimate_randoms() can be used to get a relatively noiseless estimate
        of the random coincidences.

    Currently, the randoms are estimated from the delayed coincidences using
    the following strategy:
       1. singles (one per detector) are estimated using a Maximum Likelihood
          estimator
       2. randoms-from-singles are computed per detector-pair via the usual
          product formula. These are then added together for all detector pairs
          in a certain histogram-bin in the data (accommodating for view
          mashing and axial compression).

    The actual algorithm is described in

    D. Hogg, K. Thielemans, S. Mustafovic, and T. J. Spinks,
    "A study of bias for various iterative reconstruction methods in PET,"
    in 2002 IEEE Nuclear Science Symposium Conference Record, vol. 3. IEEE,
    Nov. 2002, pp. 1519-1523 (http://dx.doi.org/10.1109/nssmic.2002.1239610).
    """

    def __init__(self, file=None):
        """init."""
        self.handle = None
        self.name = 'ListmodeToSinograms'
        if file is None:
            self.handle = pystir.cSTIR_newObject(self.name)
        else:
            self.handle = pystir.cSTIR_objectFromFile(self.name, file)
        self.output = None

    def __del__(self):
        """del."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_input(self, lm_data):
        """Sets the listmode file name, or ListmodeData object."""
        if isinstance(lm_data, str):
            parms.set_char_par(self.handle, self.name, 'input_file', lm_data)
        else:
            parms.set_parameter(self.handle, self.name, 'input', lm_data.handle)

    def set_output_prefix(self, sino_file):
        """Sets the sinograms file names prefix."""
        parms.set_char_par(self.handle, self.name, 'output', sino_file)

    def set_template(self, templ):
        """Sets the sinograms template.

        templ: either file name or AcquisitionData
        """
        if type(templ) == type('a'):
            parms.set_char_par(self.handle, self.name, 'template_file', templ)
        else:
            parms.set_parameter(self.handle, self.name, 'template', templ.handle)

    def set_time_interval(self, start, stop):
        """Sets the time interval.

        Only data scanned during this time interval will be converted.
        """
        interval = numpy.ndarray((2,), dtype=numpy.float32)
        interval[0] = start
        interval[1] = stop
        try_calling(pystir.cSTIR_setListmodeToSinogramsInterval(
            self.handle, interval.ctypes.data))

    def flag_on(self, flag):
        """Switches on (sets to 'true') a conversion flag.

        (see conversion flags description above).
        """
        try_calling(pystir.cSTIR_setListmodeToSinogramsFlag(
            self.handle, flag, 1))

    def flag_off(self, flag):
        """Switches off (sets to 'false') a conversion flag.

        (see conversion flags description above).
        """
        try_calling(pystir.cSTIR_setListmodeToSinogramsFlag(
            self.handle, flag, 0))

    def set_up(self):
        """Sets up the conversion."""
        try_calling(
            pystir.cSTIR_setupListmodeToSinogramsConverter(self.handle))

    def process(self):
        """Performs the conversion."""
        self.output = AcquisitionData()
        self.output.handle = pystir.cSTIR_convertListmodeToSinograms(
            self.handle)
        check_status(self.output.handle)

    def get_output(self):
        """Returns the sinograms as an AcquisitionData object."""
        if self.output is None:
            raise error('Conversion to sinograms not done')
        return self.output

    def estimate_randoms(self):
        """Returns an estimate of the randoms as an AcquisitionData object."""
        randoms = AcquisitionData()
        randoms.handle = pystir.cSTIR_computeRandoms(self.handle)
        check_status(randoms.handle)
        return randoms

    def get_time_at_which_num_prompts_exceeds_threshold(self, threshold):
        """Returns the time at which the number of prompts exceeds <threshold>.

        Returns -1 if no corresponding time is found.
        """
        h = pystir.cSTIR_lm_num_prompts_exceeds_threshold(
            self.handle, float(threshold))
        check_status(h, inspect.stack()[1])
        v = pyiutil.floatDataFromHandle(h)
        pyiutil.deleteDataHandle(h)
        return v


class PoissonNoiseGenerator(object):
    """
    Generates noise realisations according to Poisson statistics but allowing for scaling.

    A scaling_factor is used to multiply the input data before generating
    the Poisson random number. This means that a scaling_factor larger than 1
    will result in data with lower relative noise.

    If preserve_mean=false, the mean of the output data will
    be equal to scaling_factor*mean_of_input, otherwise it
    will be equal to mean_of_input, but then the output is no longer Poisson
    distributed.
    """

    def __init__(self, scaling_factor=1.0, preserve_mean=False):
        self.name = "PoissonNoiseGenerator"
        self.handle = pystir.cSTIR_createPoissonNoiseGenerator(scaling_factor, preserve_mean)
        check_status(self.handle)
        self.output_handle = None

    def set_seed(self, s):
        parms.set_int_par(self.handle, self.name, 'seed', s)

    def process(self, acq_data):
        self.output_handle = pystir.cSTIR_generatePoissonNoise(self.handle, acq_data.handle)
        check_status(self.output_handle)

    def get_output(self):
        if self.output_handle is None:
            raise error('Noise generating not done')
        output = AcquisitionData()
        output.handle = self.output_handle
        return output

    def generate_noisy_data(self, acq_data):
        noisy_data = AcquisitionData()
        noisy_data.handle = pystir.cSTIR_generatePoissonNoise(self.handle, acq_data.handle)
        check_status(noisy_data.handle)
        return noisy_data


class AcquisitionSensitivityModel(object):
    """
    Class that handles PET scanner detector efficiencies and attenuation.

    Is used by AcquisitionModel (see below).
    """

    def __init__(self, src, other_src=None):
        """
        Create new AcquisitionSensitivityModel object.

        Sources:
        - from a manufacturer normalisation file (supported by STIR) or
        - from ImageData object containing attenuation image (units: 1/cm) or
        - from AcquisitionData object containing bin efficiencies or
        - by chaining two existing AcquisitionSensitivityModel objects
        src: file name or ImageData object or AcquisitionData object
        other_src: AcquisitionSensitivityModel object (optional)
        """
        self.handle = None
        self.name = 'AcquisitionSensitivityModel'
        if src is None:
            return
        if isinstance(src, str):
            # create from ECAT8/GE norm file
            print('Reading manufacturer PET normalisation file from ' + src)
            handle = pyiutil.charDataHandle(src)
            self.handle = pystir.cSTIR_createPETAcquisitionSensitivityModel(
                handle, 'n')
            pyiutil.deleteDataHandle(handle)
        elif isinstance(src, ImageData):
            # create from attenuation image
            if src.handle is None:
                raise AssertionError()
            if other_src is None:
                raise AssertionError('AcquisitionSensitivityModel constructor' +
                ' with attenuation image needs an AcquisitionModel' +
                ' as second argument (for ray tracing)')
            assert_validity(other_src, AcquisitionModel)
            self.handle = pystir.cSTIR_createPETAttenuationModel(
                src.handle, other_src.handle)
        elif isinstance(src, AcquisitionData):
            # create from bin efficiencies (cf. AcquisitionModel)
            if src.handle is None:
                raise AssertionError()
            self.handle = pystir.cSTIR_createPETAcquisitionSensitivityModel(
                src.handle, 's')
        elif isinstance(src, AcquisitionSensitivityModel) and \
                isinstance(other_src, AcquisitionSensitivityModel):
            # chain two acquisition sensitivity models
            if src.handle is None:
                raise AssertionError()
            if other_src.handle is None:
                raise AssertionError()
            self.handle = pystir.cSTIR_chainPETAcquisitionSensitivityModels(
                src.handle, other_src.handle)
        else:
            raise error(
                'Wrong source in AcquisitionSensitivityModel constructor')
        check_status(self.handle)

    def set_up(self, ad):
        """Sets up the object."""
        if self.handle is None:
            raise AssertionError()
        assert_validity(ad, AcquisitionData)
        try_calling(pystir.cSTIR_setupAcquisitionSensitivityModel(
            self.handle, ad.handle))

    def normalise(self, ad):
        """Multiplies ad by the inverse n of S from AcquisitionModel (F).

        If self is a chain of two AcquisitionSensitivityModels, then n is
        a product of two normalisations.
        """
        if self.handle is None:
            raise AssertionError()
        assert_validity(ad, AcquisitionData)
        if ad.read_only:
            raise error(
                'Cannot normalise a read-only object, consider using method invert instead')
        try_calling(pystir.cSTIR_applyAcquisitionSensitivityModel(
            self.handle, ad.handle, 'normalise'))

    def unnormalise(self, ad):
        """Multiply the argument by S from AcquisitionModel (F).

        If self is a chain of two AcquisitionSensitivityModels, then S is
        a product of two un-normalisations.
        """
        if self.handle is None:
            raise AssertionError()
        assert_validity(ad, AcquisitionData)
        if ad.read_only:
            raise error(
                'Cannot unnormalise a read-only object, consider using method forward instead')
        try_calling(pystir.cSTIR_applyAcquisitionSensitivityModel(
            self.handle, ad.handle, 'unnormalise'))

    def forward(self, ad):
        """Alias of unnormalise except the argument remains unchanged

        and a new AcquisitionData equal to the argument multiplied
        by S is returned.
        """
        if self.handle is None:
            raise AssertionError()
        assert_validity(ad, AcquisitionData)
        fd = AcquisitionData()
        fd.handle = pystir.cSTIR_applyAcquisitionSensitivityModel(
            self.handle, ad.handle, 'fwd')
        check_status(fd.handle)
        return fd

    def invert(self, ad):
        """Alias of normalise except that the argument remains unchanged

        and a new AcquisitionData equal to the argument multiplied
        by the inverse of S is returned.
        """
        if self.handle is None:
            raise AssertionError()
        assert_validity(ad, AcquisitionData)
        fd = AcquisitionData()
        fd.handle = pystir.cSTIR_applyAcquisitionSensitivityModel(
            self.handle, ad.handle, 'inv')
        check_status(fd.handle)
        return fd

    @staticmethod
    def compute_attenuation_factors(sinograms, mu_map):
        '''Creates attenuation model and returns the attenuation factor (af)
        and the attenuation correction factor (acf) as AcquisitionData objects
        '''
        am = AcquisitionModelUsingRayTracingMatrix()
        attn = AcquisitionSensitivityModel(mu_map, am)
        af = AcquisitionData(sinograms)
        acf = AcquisitionData(sinograms)
        am.set_up(sinograms, mu_map)
        attn.set_up(sinograms)
        try_calling(pystir.cSTIR_computeACF(sinograms.handle, attn.handle, af.handle, acf.handle))
        return af, acf

    def __del__(self):
        """del."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)


class AcquisitionModel(object):
    """Class for PET acquisition model.

    Class for a PET acquisition model that relates an image x to the
    acquisition data y as

    (F)    y = [S] (G [P] x + [a]) + [b]

    where:
    square brackets embrace optional terms;
    P is an ImageDataProcessor, assumed to be the identity
    if not present;
    G is the geometric (ray tracing) projector from the image voxels
    to the scanner's pairs of detectors (bins);
    a and b are the additive and background terms representing
    the effects of accidental coincidences and scattering, assumed
    to be 0 if not present;
    S is the acquisition sensitivity map representing the
    detector (bin) efficiencies (see AcquisitionSensitivityModel),
    assumed to be the identity if not present.

    The computation of y for a given x by the above formula (F) is
    referred to as the forward projection, and the computation of

    (B)    z = P' G' S y

    where P' and G' are the transposes of P and G, is referred to as
    the backprojection.
    """

    def __init__(self):
        """init."""
        self.handle = None
        self.acq_templ = None
        self.img_templ = None
        # reference to the background term
        self.bt = None
        # reference to the additive term
        self.at = None
        # reference to the acquisition sensitivity model
        self.asm = None
        # default values of subset number and selected subset
        self._num_subsets = 1
        self._subset_num = 0
        # constness flag for const reference
        self.const = False

    def set_up(self, acq_templ, img_templ):
        """Sets up.

        Prepares this object for performing forward and backward
        projections;
        acq_templ:  an AcquisitionData object used as a template for
                    creating an AcquisitionData object to store forward
                    projection;
        img_templ:  an ImageData object used as a template for creating an
                    ImageData object to store backward projection.
        """
        if self.const:
            raise RuntimeError('cannot set_up const object')

        assert_validity(acq_templ, AcquisitionData)
        assert_validity(img_templ, ImageData)

        try_calling(pystir.cSTIR_setupAcquisitionModel(
            self.handle, acq_templ.handle, img_templ.handle))

        self.acq_templ = acq_templ
        self.img_templ = img_templ

    def norm(self, subset_num=0, num_subsets=1, num_iter=2, verb=0):
        """Computes the norm of a part the linear operator S G.

        See the docstring for method forward for the desciption of
        the partial forward projection and the arguments.
        """
        assert self.handle is not None
        handle = pystir.cSTIR_acquisitionModelNorm \
                 (self.handle, subset_num, num_subsets, num_iter, verb)
        check_status(handle)
        r = pyiutil.floatDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return r;

    def set_additive_term(self, at):
        """Sets additive term.

        Sets the additive term a in the acquisition model;
        at:  an AcquisitionData object containing a.
        """
        if self.const:
            raise RuntimeError('cannot set_additive_term for a const object')
        assert_validity(at, AcquisitionData)
        parms.set_parameter(
            self.handle, 'AcquisitionModel', 'additive_term', at.handle)
        # save reference to the additive term
        self.at = at

    def set_background_term(self, bt):
        """
        Sets the background term b in the acquisition model.

        bt:  an AcquisitionData object containing b.
        """
        if self.const:
            raise RuntimeError('cannot set_background_term for a const object')
        assert_validity(bt, AcquisitionData)
        parms.set_parameter(
            self.handle, 'AcquisitionModel', 'background_term', bt.handle)
        # save reference to the background term
        self.bt = bt

    def set_image_data_processor(self, processor):
        """
        Sets the image processor P in the acquisition model.

        processor:  an ImageDataProcessor object.
        """
        if self.const:
            msg = 'cannot set_image_data_processor'
            raise RuntimeError(msg + ' for a const object')
        assert_validity(processor, ImageDataProcessor)
        parms.set_parameter(
            self.handle, 'AcquisitionModel', 'image_data_processor',
            processor.handle)

    def get_background_term(self):
        """Returns the background term b of the AcquisitionModel (F).
        """
        if self.bt is None:
            if self.acq_templ is None:
                raise RuntimeError('AcquisitionModel.set_up() call missing')
            self.bt = AcquisitionData(self.acq_templ)
            self.bt.fill(0)
        return self.bt

    def get_additive_term(self):
        """Returns the additive term a of the AcquisitionModel (F).
        """
        if self.at is None:
            if self.acq_templ is None:
                raise RuntimeError('AcquisitionModel.set_up() call missing')
            self.at = AcquisitionData(self.acq_templ)
            self.at.fill(0)
        return self.at

    def get_constant_term(self):
        """Returns the constant term S a + b in (F).
        """
        if self.asm is not None:
            return self.asm.forward(self.get_additive_term()) + \
                   self.get_background_term()
        else:
            return self.get_additive_term() + self.get_background_term()

    def set_acquisition_sensitivity(self, asm):
        """Sets S in the AcquisitionModel (F).

        asm:  an AcquisitionSensitivityModel object.
        """
        if self.const:
            msg = 'cannot set_acquisition_sensitivity'
            raise RuntimeError(msg + ' for a const object')
        assert_validity(asm, AcquisitionSensitivityModel)
        parms.set_parameter(
            self.handle, 'AcquisitionModel', 'asm', asm.handle)
        # save reference to the Acquisition Sensitivity Model
        self.asm = asm

    def forward(self, image, subset_num=None, num_subsets=None, out=None):
        """Returns the [partial] forward projection of image.

        If <num_subsets> is greater than 1, then the forward projection data
        y in AcquisitionModel (F) is divided into <num_subsets> subsets and
        only subset <subset_num> is computed; otherwise the full forward
        projection is computed.
        image      : an ImageData object.
        subset_num : int, optional
                     subset number to forward project to; if None, subset_num
                     is set to self.subset_num.
        num_subsets: int, optional
                     number of subsets y is divided into; if None, num_subsets
                     is set to self.num_subsets.
        out        : an existing AcquisitionData object, optional
                     the destination for the projection; if None a new
                     AcquisitionData object will be returned.
        """
        assert_validity(image, ImageData)
        if subset_num is None:
            subset_num = self.subset_num
        if num_subsets is None:
            num_subsets = self.num_subsets
        if out is None:
            ad = AcquisitionData()
            ad.handle = pystir.cSTIR_acquisitionModelFwd(
                self.handle, image.handle, subset_num, num_subsets)
            check_status(ad.handle)
            return ad
        ad = out
        assert_validity(ad, AcquisitionData)
        try_calling(pystir.cSTIR_acquisitionModelFwdReplace(
            self.handle, image.handle, subset_num, num_subsets, ad.handle))

    def backward(self, ad, subset_num=None, num_subsets=None, out=None):
        """
        Return the [partial] backward projection of ad.

        ad         : an AcquisitionData object.
        subset_num : int, optional
                     subset number to backproject; if None, is set to
                     self.subset_num.
        num_subsets: int, optional
                     the number of subsets of ad; if None, is set to
                     self.num_subsets.
        out        : optional AcquisitionData to store the result into.
                     Default None, if None a new AcquisitionData will be
                     returned.
        """
        assert_validity(ad, AcquisitionData)
        if subset_num is None:
            subset_num = self.subset_num
        if num_subsets is None:
            num_subsets = self.num_subsets
        if out is None:
            image = ImageData()
            image.handle = pystir.cSTIR_acquisitionModelBwd(
                self.handle, ad.handle, subset_num, num_subsets)
            check_status(image.handle)
            return image
        assert_validity(out, ImageData)
        try_calling(pystir.cSTIR_acquisitionModelBwdReplace(
                self.handle, ad.handle, subset_num, num_subsets, out.handle))

    def get_linear_acquisition_model(self):
        """Returns the linear part L = S G P of self.
        """
        am = AcquisitionModel()
        am.handle = pystir.cSTIR_linearAcquisitionModel(self.handle)
        check_status(am.handle)
        am.const = True # am to be a const reference of self
        return am

    def direct(self, image, out=None):
        '''Projects an image into the (simulated) acquisition space,
           calls forward with num_subset and on subset_num members

           Added for CCPi CIL compatibility
           https://github.com/CCPPETMR/SIRF/pull/237#issuecomment-439894266
        '''
        if self.is_linear():
            return self.forward(image, \
                                subset_num=self.subset_num, \
                                num_subsets=self.num_subsets, \
                                out=out)
        else:
            raise error('AcquisitionModel is not linear\nYou can get the ' +
                        'linear part of the AcquisitionModel with ' +
                        'get_linear_acquisition_model')

    def adjoint(self, ad, out=None):
        '''Back-projects acquisition data into image space, if the
           AcquisitionModel is linear

           calls backward with num_subset and on subset_num members
           Added for CCPi CIL compatibility
           https://github.com/CCPPETMR/SIRF/pull/237#issuecomment-439894266
        '''
        if self.is_linear():
            return self.backward(ad, subset_num=self.subset_num,
                                 num_subsets=self.num_subsets, out=out)
        else:
            raise error('AcquisitionModel is not linear. You can get the '
                        'linear part of the AcquisitionModel with '
                        'get_linear_acquisition_model')

    def is_affine(self):
        """Returns True if the acquisition model is affine.
        """
        return True

    def is_linear(self):
        """Return true if the acquisition model constant term is zero.

        i.e. corresponding to
        A*x, with zero background term.
        """
        if self.bt is None and self.at is None:
            return True
        else:
            if self.bt is None and self.at is not None:
                return self.at.norm() == 0
            elif self.bt is not None and self.at is None:
                return self.bt.norm() == 0
            else:
                return self.bt.norm() == 0 and self.at.norm() == 0

    def range_geometry(self):
        """Return the template of AcquisitionData."""
        geom = AcquisitionData()
        geom.handle = pystir.cSTIR_parameter(
            self.handle, 'AcquisitionModel', 'range geometry')
        check_status(geom.handle)
        return geom

    def domain_geometry(self):
        """Return the template of ImageData."""
        geom = ImageData()
        geom.handle = pystir.cSTIR_parameter(
            self.handle, 'AcquisitionModel', 'domain geometry')
        check_status(geom.handle)
        return geom

    @property
    def subset_num(self):
        '''Selected subset number

        This value is used by direct and adjoint methods and are the
        default values used by forward and back projection for their
        parameter subset_num.

        Default value is 0.
        '''
        return self._subset_num

    @property
    def num_subsets(self):
        '''Number of subsets to divide the AcquisitionData during projection

        This value is used by the direct and adjoint methods. Additionally,
        this value is the default value used by forward and back projection
        for the parameter num_subsets.

        Default value is 1 and corresponds to forward/backward projecting
        the whole dataset.

        '''
        return self._num_subsets

    @subset_num.setter
    def subset_num(self, value):
        '''setter for subset_num

        value: int >= 0 and < num_subsets
        '''
        if isinstance (value, Integral):
            if value < self.num_subsets and value >= 0:
                self._subset_num = value
            else:
                raise ValueError("Expected a subset number below {} and larger or equal than 0. Got {}"\
                    .format(self.subset_num, value))
        else:
            raise ValueError("Expected an integer. Got {}".format(type(value)))

    @num_subsets.setter
    def num_subsets(self, value):
        '''setter for num_subsets

        value: int > 0.
        Allows to set the number of subsets the AcquisitionModel operates on.
        Notice that reassigning the num_subsets to any valid number will also
        set the property subset_num to 0.
        '''
        if isinstance (value, Integral):
            if value > 0:
                self._num_subsets = value
                self.subset_num = 0

            else:
                raise ValueError("Expected a subset number larger than 0. Got {}"\
                    .format(value))
        else:
            raise ValueError("Expected an integer. Got {}".format(type(value)))



class AcquisitionModelUsingMatrix(AcquisitionModel):
    """PET acquisition model with sparse matrix.

    Class for a PET acquisition model that uses (implicitly) a sparse
    matrix for G in AcquisitionModel (F).
    """

    def __init__(self, matrix=None):
        """Creates an AcquisitionModelUsingMatrix object,

        optionally setting the ray tracing matrix to be used for projecting.
        matrix: an object to represent G in (F).
        # TODO will need to try different matrices here
        """
        super(AcquisitionModelUsingMatrix, self).__init__()
        self.handle = None
        self.name = 'AcqModUsingMatrix'
        self.handle = pystir.cSTIR_newObject(self.name)
        check_status(self.handle)
        if matrix is None:
            matrix = RayTracingMatrix()
        self.set_matrix(matrix)

    def __del__(self):
        """del."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_matrix(self, matrix):
        '''
        Sets the matrix G to be used for projecting;
        matrix:  a matrix object to represent G in acquisition model (F).
        '''
        # The following allows for different matrices
        try:
            parms.set_parameter(self.handle, self.name, 'matrix', matrix.handle)
        except:
            raise AssertionError('Unknown matrix type.')


class AcquisitionModelUsingRayTracingMatrix(AcquisitionModelUsingMatrix):
    """PET acquisition model with RayTracingMatrix.

    Class for a PET acquisition model that uses (implicitly) a RayTracingMatrix
    for G in AcquisitionModel (F).
    """

    def __init__(self, matrix=None):
        """Create an AcquisitionModelUsingMatrix object,

        optionally setting the ray tracing matrix to be used for projecting;
        matrix: a RayTracingMatrix object to represent G in (F).
        """
        self.handle = None
        if matrix is None:
            matrix = RayTracingMatrix()
        assert_validity(matrix, RayTracingMatrix)
        super(AcquisitionModelUsingRayTracingMatrix, self).__init__(matrix)

    def __del__(self):
        """del."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_matrix(self, matrix):
        """Sets the ray tracing matrix to be used for projecting.

        matrix: a RayTracingMatrix object to represent G in (F).
        """
        assert_validity(matrix, RayTracingMatrix)
        parms.set_parameter(self.handle, self.name, 'matrix', matrix.handle)

    def get_matrix(self):
        """Returns the ray tracing matrix used for projecting.

        matrix: a RayTracingMatrix object representing G in (F).
        """
        matrix = RayTracingMatrix()
        matrix.handle = pystir.cSTIR_parameter(
            self.handle, self.name, 'matrix')
        check_status(matrix.handle)
        return matrix

    def set_num_tangential_LORs(self, value):
        """See :func:`~sirf.STIR.RayTracingMatrix.set_num_tangential_LORs`."""
        return self.get_matrix().set_num_tangential_LORs(value)

    def get_num_tangential_LORs(self):
        """See :func:`~sirf.STIR.RayTracingMatrix.get_num_tangential_LORs`."""
        return self.get_matrix().get_num_tangential_LORs()


if SIRF_HAS_NiftyPET:
    class AcquisitionModelUsingNiftyPET(AcquisitionModel):
        """PET acquisition model that uses NiftyPET projector.

        Class for a PET acquisition model that uses the NiftyPET projector
        for G in AcquisitionModel (F).
        """

        def __init__(self):
            """Create an AcquisitionModelUsingNiftyPET object."""
            super(AcquisitionModelUsingNiftyPET, self).__init__()
            self.name = 'AcqModUsingNiftyPET'
            self.handle = pystir.cSTIR_newObject(self.name)
            check_status(self.handle)

        def __del__(self):
            """del."""
            if self.handle is not None:
                pyiutil.deleteDataHandle(self.handle)

        def set_cuda_verbosity(self, verbosity):
            """Set the verbosity of the CUDA code."""
            if verbosity:
                v = 1
            else:
                v = 0
            parms.set_int_par(self.handle, self.name, 'cuda_verbosity', v)

        def set_use_truncation(self, use_truncation):
            """Set use truncation.

            Whether or not to truncate FOV before forward- and
            after back-projection.
            """
            if use_truncation:
                v = 1
            else:
                v = 0
            parms.set_int_par(self.handle, self.name, 'use_truncation', v)


if SIRF_HAS_Parallelproj:
    class AcquisitionModelUsingParallelproj(AcquisitionModel):
        """PET acquisition model that uses Parallelproj projector.

        Class for a PET acquisition model that uses the Parallelproj projector
        (see https://github.com/gschramm/parallelproj) for G in
        AcquisitionModel (F).
        """

        def __init__(self):
            """Create an AcquisitionModelUsingNiftyPET object."""
            super(AcquisitionModelUsingParallelproj, self).__init__()
            self.name = 'AcqModUsingParallelproj'
            self.handle = pystir.cSTIR_newObject(self.name)
            check_status(self.handle)

        def __del__(self):
            """del."""
            if self.handle is not None:
                pyiutil.deleteDataHandle(self.handle)


class Prior(object):
    """Class for objects handling the prior, a penalty term

    to be added to the objective function maximized by iterative
    reconstruction algorithms.
    """

    def __init__(self):
        """init."""
        self.handle = None

    def __del__(self):
        """del."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def __call__(self, image):
        '''Returns the prior value on the specified image (alias of value()).

        image: ImageData object'''
        return self.value(image)

    def set_penalisation_factor(self, value):
        """Sets penalisation factor.

        Sets the factor by which the penalty term (prior) is to be multiplied
        before adding to the objective function.
        """
        parms.set_float_par(
            self.handle, 'GeneralisedPrior', 'penalisation_factor', value)
        return self

    def get_penalisation_factor(self):
        """Returns the penalty factor in front of the prior."""
        return parms.float_par(
            self.handle, 'GeneralisedPrior', 'penalisation_factor')

    def get_value(self, image):
        """Returns the value of the prior.

        Returns the value of the prior for the specified image.
        image: ImageData object
        """
        assert_validity(image, ImageData)
        handle = pystir.cSTIR_priorValue(self.handle, image.handle)
        check_status(handle)
        v = pyiutil.doubleDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return v

    def value(self, image):
        """Returns the value of the prior (alias of get_value())."""
        return self.get_value(image)

    def get_gradient(self, image, out=None, **kwargs):
        """Returns gradient of the prior.

        Returns the value of the gradient of the prior for the specified image.
        image: ImageData object
        """
        assert_validity(image, ImageData)
        if out is None:
            out = ImageData()
        if out.handle is None:
            out.handle = pystir.cSTIR_priorGradient(self.handle, image.handle)
        else:
            assert_validities(image, out)
            pystir.cSTIR_computePriorGradient(self.handle, image.handle, out.handle)
        check_status(out.handle)
        return out

    def gradient(self, image, out=None, **kwargs):
        """Returns the gradient of the prior (alias of get_gradient())."""

        return self.get_gradient(image, out)

    def set_up(self, image):
        """Sets up."""
        try_calling(pystir.cSTIR_setupPrior(self.handle, image.handle))

    def accumulate_Hessian_times_input(self, current_estimate, input_, out=None):
        """Computes the multiplication of the Hessian with a vector and adds it to output.
        """
        if out is None or out.handle is None:
            out = input_.get_uniform_copy(0.0)
        try_calling(pystir.cSTIR_priorAccumulateHessianTimesInput
            (self.handle, current_estimate.handle, input_.handle, out.handle))
        return out

    def multiply_with_Hessian(self, current_estimate, input_, out=None):
        """Computes the multiplication of the Hessian at current_estimate with a vector.
        """
        if out is None or out.handle is None:
            out = input_.get_uniform_copy(0.0)
        try_calling(pystir.cSTIR_priorComputeHessianTimesInput
            (self.handle, current_estimate.handle, input_.handle, out.handle))
        return out

class QuadraticPrior(Prior):
    r"""Class for the prior that is a quadratic function of the image values.

    Implements a quadratic Gibbs prior.

    The gradient of the prior for the image lambda is computed at voxel r as
    the sum of

        delta(r, s) = w(s - r)*kappa(r)*kappa(s)*(lambda(r) - lambda(s))

    over all voxels s where the weight w(s - r) is non-zero, kappa being an
    image used to have spatially-varying penalties such as in Jeff Fessler's
    papers.It should have identical dimensions to the image for which the
    penalty is computed. If kappa is not set, this class will
    effectively use 1 for all kappa's.

    By default, a 3x3 or 3x3x3 neigbourhood is used where the weights are set
    to x-voxel_size divided by the Euclidean distance between the points.
    """

    def __init__(self):
        """init."""
        self.handle = None
        self.name = 'QuadraticPrior'
        self.handle = pystir.cSTIR_newObject(self.name)
        check_status(self.handle)

    def __del__(self):
        """del."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_kappa(self, image):
        """Sets kappa."""
        assert_validity(image, ImageData)
        parms.set_parameter(self.handle, 'QuadraticPrior', 'kappa', image.handle)

    def get_kappa(self):
        """Returns kappa."""
        image = ImageData()
        image.handle = pystir.cSTIR_parameter(self.handle, 'QuadraticPrior', 'kappa')
        check_status(image.handle)
        return image


class LogcoshPrior(Prior):
    r"""Class for Log-cosh Prior.

    Implements the prior, Log-cosh Prior, one of the earliest uses in P. J.
    Green's paper "Bayesian reconstructions from emission tomography data using
    a modified EM algorithm," in IEEE Transactions on Medical Imaging, vol. 9,
    no. 1, pp. 84-93, March 1990, doi: 10.1109/42.52985.

    The prior has one parameter the scalar, it is the edge-preservation parameter.

    The log-cosh function is given by:
    \f[
        f = \sum_{r,dr} w_{dr} \frac{1}{2 s^2}  log(cosh(s(\lambda_r - \lambda_{r+dr}))) * \kappa_r * \kappa_{r+dr}$
    \f]

    Kappa is a spatially varying penalty strength.
    """

    def __init__(self):
        """init."""
        self.handle = None
        self.name = 'LogcoshPrior'
        self.handle = pystir.cSTIR_newObject(self.name)
        check_status(self.handle)

    def __del__(self):
        """del."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_scalar(self, v):
        """Sets scalar."""
        parms.set_float_par(self.handle, 'LogcoshPrior', 'scalar', v)

    def get_scalar(self):
        """Returns scalar."""
        return parms.float_par(self.handle, 'LogcoshPrior', 'scalar')

    def set_kappa(self, image):
        """Sets kappa."""
        assert_validity(image, ImageData)
        parms.set_parameter(self.handle, 'LogcoshPrior', 'kappa', image.handle)

    def get_kappa(self):
        """Returns kappa."""
        image = ImageData()
        image.handle = pystir.cSTIR_parameter(self.handle, 'LogcoshPrior', 'kappa')
        check_status(image.handle)
        return image


class RelativeDifferencePrior(Prior):
    r"""Class for Relative Difference Prior.

    Implements the prior, Relative Difference Prior, proposed by Johan Nuyts et.
    al in "A concave prior penalizing relative differences for
    maximum-a-posteriori reconstruction in emission tomography," in IEEE
    Transactions on Nuclear Science, vol. 49, no. 1, pp. 56-60, Feb. 2002,
    doi: 10.1109/TNS.2002.998681.

    The value of the prior is computed as follows:

    \f[
    f = \sum_{r,dr} \frac{w_{dr}}{2} \frac{(\lambda_r - \lambda_{r+dr})^2}{(\lambda_r+ \lambda_{r+dr} + \gamma |\lambda_r - \lambda_{r+dr}| + \epsilon)} * \kappa_r * \kappa_{r+dr}
    \f]

    The prior has 2 parameters epsilon and gamma. The former is to ensure
    numerical stability and the gamma is the edge-preservation parameters
    typically set as 2 in clinical practice (citation required).

    Kappa is a spatially varying penalty strength.
    """

    def __init__(self):
        """init."""
        self.handle = None
        self.name = 'RelativeDifferencePrior'
        self.handle = pystir.cSTIR_newObject(self.name)
        check_status(self.handle)

    def __del__(self):
        """del."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_gamma(self, v):
        """Sets gamma."""
        parms.set_float_par(self.handle, 'RelativeDifferencePrior', 'gamma', v)

    def get_gamma(self):
        """Returns gamma."""
        return parms.float_par(self.handle, 'RelativeDifferencePrior', 'gamma')

    def set_epsilon(self, v):
        """Sets epsilon."""
        parms.set_float_par(self.handle, 'RelativeDifferencePrior', 'epsilon', v)

    def get_epsilon(self):
        """Returns epsilon."""
        return parms.float_par(self.handle, 'RelativeDifferencePrior', 'epsilon')

    def set_kappa(self, image):
        """Sets kappa."""
        assert_validity(image, ImageData)
        parms.set_parameter(self.handle, 'RelativeDifferencePrior', 'kappa', image.handle)

    def get_kappa(self):
        """Returns kappa."""
        image = ImageData()
        image.handle = pystir.cSTIR_parameter(self.handle, 'RelativeDifferencePrior', 'kappa')
        check_status(image.handle)
        return image

if STIR_WITH_CUDA:
    class CudaRelativeDifferencePrior(RelativeDifferencePrior):
        r"""Class for Relative Difference Prior using CUDA computations

        Identical to RelativeDifferencePrior, but using STIR's CUDA implementation.
        """

        def __init__(self):
            """init."""
            self.name = 'CudaRelativeDifferencePrior'
            self.handle = pystir.cSTIR_newObject(self.name)
            check_status(self.handle)

        def __del__(self):
            """del."""
            if self.handle is not None:
                pyiutil.deleteDataHandle(self.handle)


class PLSPrior(Prior):
    r"""Class for Parallel Level Sets prior.

    Implements the anatomical penalty function, Parallel Level Sets (PLS),
    proposed by Matthias J. Ehrhardt et. al in "PET Reconstruction With an
    Anatomical MRI Prior Using Parallel Level Sets", IEEE Trans. med. Imag.,
    vol. 35, no. 9, Sep 2016 (https://doi.org/10.1109/TMI.2016.2549601).
    Note that PLS becomes smoothed TV when a uniform anatomical image is
    provided.

    The prior has 2 parameters alpha and eta. It is computed for an image \f$
    f f$ as

        phi(f) = sqrt(alpha^2 + |grad f|^2 - <grad f, xi>^2)

    where f is the PET image, alpha controls the edge-preservation property
    of PLS, and depends on the scale of the emission image, xi is the
    normalised gradient of the anatomical image calculated as follows:

        xi = 1/sqrt(|grad v|^2 + eta^2) grad v

    where v is the anatomical image, and eta safeguards against the division
    by zero and depends on the scale of the anatomical image.

    """

    def __init__(self):
        """init."""
        self.handle = None
        self.name = 'PLSPrior'
        self.handle = pystir.cSTIR_newObject(self.name)
        check_status(self.handle)

    def __del__(self):
        """del."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_only_2D(self, tf):
        """Sets only_2D parameter."""
        v = 1 if tf else 0
        parms.set_int_par(self.handle, 'PLSPrior', 'only_2D', v)

    def get_only_2D(self):
        """Returns only_2D value."""
        v = parms.int_par(self.handle, 'PLSPrior', 'only_2D')
        return v != 0

    def set_alpha(self, v):
        """Sets alpha."""
        parms.set_float_par(self.handle, 'PLSPrior', 'alpha', v)

    def get_alpha(self):
        """Returns alpha."""
        return parms.float_par(self.handle, 'PLSPrior', 'alpha')

    def set_eta(self, v):
        """Sets eta."""
        parms.set_float_par(self.handle, 'PLSPrior', 'eta', v)

    def get_eta(self):
        """Returns eta."""
        return parms.float_par(self.handle, 'PLSPrior', 'eta')

    def set_anatomical_image(self, image):
        """Sets anatomical image."""
        assert_validity(image, ImageData)
        parms.set_parameter(
            self.handle, 'PLSPrior', 'anatomical_image', image.handle)

    def get_anatomical_image(self):
        """Returns anatomical image."""
        image = ImageData()
        image.handle = pystir.cSTIR_parameter(
            self.handle, 'PLSPrior', 'anatomical_image')
        check_status(image.handle)
        return image

    def get_anatomical_grad(self, direction, out=None):
        """Returns anatomical gradient."""
        if out is None:
            image = ImageData()
        else:
            image = out
        image.handle = pystir.cSTIR_PLSPriorAnatomicalGradient(self.handle, direction)
        check_status(image.handle)
        return image

    def set_anatomical_filename(self, filename):
        """Sets anatomical filename."""
        parms.set_char_par(
            self.handle, 'PLSPrior', 'anatomical_filename', filename)

    def set_kappa(self, image):
        """Sets kappa."""
        assert_validity(image, ImageData)
        parms.set_parameter(self.handle, 'PLSPrior', 'kappa', image.handle)

    def get_kappa(self):
        """Returns kappa."""
        image = ImageData()
        image.handle = pystir.cSTIR_parameter(self.handle, 'PLSPrior', 'kappa')
        check_status(image.handle)
        return image

    def set_kappa_filename(self, filename):
        """Sets kappa filename."""
        parms.set_char_par(self.handle, 'PLSPrior', 'kappa_filename', filename)

    def get_norm(self):
        """Returns norm."""
        image = ImageData()
        image.handle = pystir.cSTIR_parameter(self.handle, 'PLSPrior', 'norm')
        check_status(image.handle)
        return image


class ObjectiveFunction(object):
    """Class for the objective function

    maximised by the iterative reconstruction algorithms.
    """

    def __init__(self):
        """init."""
        self.handle = None

    def __del__(self):
        """del."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_prior(self, prior):
        """Sets the prior,

        a penalty term to be added to the objective function."""
        assert_validity(prior, Prior)
        parms.set_parameter(self.handle, 'GeneralisedObjectiveFunction',
                            'prior', prior.handle)
        self.prior = prior

    def get_prior(self):
        """Returns the prior currently used by this objective function."""
        prior = Prior()
        prior.handle = pystir.cSTIR_parameter(
            self.handle, 'GeneralisedObjectiveFunction', 'prior')
        check_status(prior.handle)
        return prior

    def set_num_subsets(self, n):
        """Sets the number of subsets.

        Sets the number of subsets of ray projections to be used for computing
        additive components of the gradient used by Ordered Subset algorithms
        for maximizing this objective function.
        If the ray tracing projector G is a matrix, the subsets in question are
        subsets of its rows.
        n: number of subsets, Python integer scalar
        """
        parms.set_int_par(
            self.handle, 'GeneralisedObjectiveFunction', 'num_subsets', n)

    def get_num_subsets(self):
        """Returns the number of subsets."""
        return parms.int_par(
            self.handle, 'GeneralisedObjectiveFunction', 'num_subsets')

    def set_up(self, image):
        """Prepares this object for use.

        image: ImageData object
        """
        assert_validity(image, ImageData)
        try_calling(pystir.cSTIR_setupObjectiveFunction(
            self.handle, image.handle))

    def value(self, image):
        """Returns the value of this objective function on the specified image.

        image: ImageData object
        """
        assert_validity(image, ImageData)
        handle = pystir.cSTIR_objectiveFunctionValue(self.handle, image.handle)
        check_status(handle)
        v = pyiutil.doubleDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return v

    def __call__(self, image):
        '''Alias of value: Returns the value of this objective function on the specified image.

        image: ImageData object'''
        return self.value(image)

    def get_value(self, image):
        """Returns the value of this objective function on the specified image.

        image: ImageData object
        """
        return self.value(image)

    def gradient(self, image, subset=-1, out=None):
        """Returns the value of the additive component of the gradient

        of this objective function on the specified image corresponding to the
        specified subset (see set_num_subsets() method).
        If no subset is specified, returns the full gradient, i.e. the sum of
        the subset components.

        Parameters:

        image: ImageData object
        subset: Python integer scalar, optional, default -1
                If subset is -1 it returns the full gradient, otherwise the
                gradient component corresponding to the specified subset.
        out: ImageData object, optional, default None
                the destination for the gradient; if None a new ImageData object
                will be returned. If 'out' is the same as 'image', the result will 
                be stored in a temporary object and then copied back to 'image', and
                returned.
        """
        assert_validity(image, ImageData)
        if out is None:
            out = ImageData()
        inline = False
        if out.handle == image.handle:
            out = ImageData()
            inline = True

        if out.handle is None:
            out.handle = pystir.cSTIR_objectiveFunctionGradient(self.handle, image.handle, subset)
        else:
            assert_validities(image, out)
            pystir.cSTIR_computeObjectiveFunctionGradient(self.handle, image.handle, subset, out.handle)
        check_status(out.handle)
        if inline:
            image.fill(out)
            return image
        return out

    def get_gradient(self, image, out=None):
        """Returns the gradient of the objective function on specified image.

        image: ImageData object
        """
        return self.gradient(image, -1, out)

    def get_subset_gradient(self, image, subset, out=None):
        """Returns the value of the additive component of the gradient

        of this objective function on <image> corresponding to the specified
        subset (see set_num_subsets() method).
        image: ImageData object
        subset: Python integer scalar
        """
        return self.gradient(image, subset, out)

    def accumulate_Hessian_times_input(self, current_estimate, input_, subset=-1, out=None):
        """Computes the multiplication of the Hessian at current_estimate with a vector and adds it to output.
        """
        if out is None or out.handle is None:
            out = input_.get_uniform_copy(0.0)
        try_calling(pystir.cSTIR_objectiveFunctionAccumulateHessianTimesInput
            (self.handle, current_estimate.handle, input_.handle, subset, out.handle))
        return out

    def multiply_with_Hessian(self, current_estimate, input_, subset=-1, out=None):
        """Computes the multiplication of the Hessian at current_estimate with a vector.
        """
        if out is None or out.handle is None:
            out = input_.get_uniform_copy(0.0) # actual value doesn't matter as STIR will overwrite it.
        try_calling(pystir.cSTIR_objectiveFunctionComputeHessianTimesInput
            (self.handle, current_estimate.handle, input_.handle, subset, out.handle))
        return out

    @abc.abstractmethod
    def get_subset_sensitivity(self, subset):
        #print('in base class ObjectiveFunction')
        pass


class PoissonLogLikelihoodWithLinearModelForMean(ObjectiveFunction):
    """Class for STIR PoissonLogLikelihoodWithLinearModelForMean object.

    See:
    http://stir.sourceforge.net/documentation/doxy/html/classstir_1_1PoissonLogLikelihoodWithLinearModelForMean.html
    """

    def __init__(self):
        """init."""
        self.handle = None

    def __del__(self):
        """del."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
#    def set_sensitivity_filename(self, name):
#        parms.set_char_par\
#            (self.handle, 'PoissonLogLikelihoodWithLinearModelForMean',\
#             'sensitivity_filename', name)
#    def set_use_subset_sensitivities(self, flag):
#        parms.set_char_par\
#            (self.handle, 'PoissonLogLikelihoodWithLinearModelForMean',\
#             'use_subset_sensitivities', repr(flag))

    def set_recompute_sensitivity(self, flag):
        """Sets recompute sensitivity flag."""
        parms.set_char_par(
            self.handle, 'PoissonLogLikelihoodWithLinearModelForMean',
            'recompute_sensitivity', repr(flag))

    def get_subset_sensitivity(self, subset):
        """Returns subset sensitivity.

        Returns an ImageData object containing sensitivity image for the
        specified subset.
        """
        ss = ImageData()
        ss.handle = pystir.cSTIR_subsetSensitivity(self.handle, subset)
        check_status(ss.handle)
        return ss

    def get_backprojection_of_acquisition_ratio(self, image, subset, out=None):
        """Returns backprojection of measured to estimated acquisition ratio.

        Returns the back-projection of the ratio of the measured and estimated
        acquisition data.
        """
        assert_validity(image, ImageData)
        if out is None:
            out = ImageData()
        if out.handle is None:
            out.handle = pystir.cSTIR_objectiveFunctionGradientNotDivided(
                self.handle, image.handle, subset)
        else:
            assert_validities(image, out)
            pystir.cSTIR_computeObjectiveFunctionGradientNotDivided(
                self.handle, image.handle, subset, out.handle)
        check_status(out.handle)
        return out


class PoissonLogLikelihoodWithLinearModelForMeanAndProjData(
        PoissonLogLikelihoodWithLinearModelForMean):
    """Class for STIR type of Poisson loglikelihood object.

    Specifically, PoissonLogLikelihoodWithLinearModelForMeanAndProjData. See:
    http://stir.sourceforge.net/documentation/doxy/html/classstir_1_1PoissonLogLikelihoodWithLinearModelForMeanAndProjData.html
    """

    def __init__(self):
        """init."""
        self.handle = None
        self.name = 'PoissonLogLikelihoodWithLinearModelForMeanAndProjData'
        self.handle = pystir.cSTIR_newObject(self.name)
        check_status(self.handle)

    def __del__(self):
        """del."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_input_filename(self, name):
        """Sets the name of the file containing acquisition data."""
        parms.set_char_par(
            self.handle, self.name, 'input_filename', name)
#    def set_zero_seg0_end_planes(self, flag):
#        parms.set_char_par\
#            (self.handle, self.name, 'zero_seg0_end_planes', repr(flag))
#    def set_max_segment_num_to_process(self, n):
#        parms.set_int_par(
#           self.handle, self.name, 'max_segment_num_to_process', n)

    def set_acquisition_model(self, am):
        """Sets the acquisition model to be used by this objective function."""
        assert_validity(am, AcquisitionModel)
        parms.set_parameter(
            self.handle, self.name, 'acquisition_model', am.handle)
#    def get_acquisition_model(self):
#        """
#        Returns the acquisition model used by this objective function.
#        """
#        am = AcquisitionModelUsingMatrix()
#        if am.handle is not None:
#            pyiutil.deleteDataHandle(am.handle)
#        am.handle = pystir.cSTIR_parameter\
#            (self.handle, self.name, 'acquisition_model')
#        check_status(am.handle)
#        return am

    def set_acquisition_data(self, ad):
        """Sets the acquisition data to be used by this objective function."""
        assert_validity(ad, AcquisitionData)
        parms.set_parameter(
            self.handle, self.name, 'acquisition_data', ad.handle)


class PoissonLogLikelihoodWithLinearModelForMeanAndListModeDataWithProjMatrixByBin(PoissonLogLikelihoodWithLinearModelForMean):
#(ObjectiveFunction):
    """Class for a STIR type of Poisson loglikelihood object for listmode data.

    Specifically, PoissonLogLikelihoodWithLinearModelForMeanAndListModeDataWithProjMatrixByBin.
    """

    def __init__(self):
        """init."""
        self.handle = None
        self.name = 'PoissonLogLikelihoodWithLinearModelForMeanAndListModeDataWithProjMatrixByBin'
        self.handle = pystir.cSTIR_newObject(self.name)
        check_status(self.handle)

    def __del__(self):
        """del."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_cache_path(self, path):
        parms.set_char_par(self.handle, self.name, 'cache_path', path)

    def get_cache_path(self):
        return parms.char_par(self.handle, self.name, 'cache_path')

    def set_time_interval(self, start, stop):
        """Sets the time interval.

        Only data scanned during this time interval will be converted.
        """
        interval = numpy.ndarray((2,), dtype=numpy.float32)
        interval[0] = start
        interval[1] = stop
        try_calling(pystir.cSTIR_objFunListModeSetInterval(
            self.handle, interval.ctypes.data))

    def set_acquisition_data(self, ad):
        assert_validity(ad, ListmodeData)
        parms.set_parameter(
            self.handle, self.name, 'acquisition_data', ad.handle)

    def set_acquisition_model(self, am):
        """Sets the acquisition model to be used by this objective function."""
        assert_validity(am, AcquisitionModel)
        parms.set_parameter(
            self.handle, self.name, 'acquisition_model', am.handle)

    # disabled for now as this doesn't work yet in STIR
    # def set_skip_lm_input_file(self, tf):
    #    flag = 1 if tf else 0
    #    parms.set_int_par(
    #            self.handle, self.name, 'skip_lm_input_file', flag)

    def set_skip_balanced_subsets(self, tf):
        ''' if tf=True, disable the check for balanced subsets (somewhat dangerous)'''
        flag = 1 if tf else 0
        parms.set_int_par(
                self.handle, self.name, 'skip_balanced_subsets', flag)

    def set_max_segment_num_to_process(self, diff):
        parms.set_int_par(
                self.handle, self.name, 'max_segment_num_to_process', diff)

    def set_recompute_cache(self, tf):
        ''' if tf=True (and cache_size>0, recompute the listmode cache'''
        flag = 1 if tf else 0
        parms.set_int_par(
                self.handle, self.name, 'recompute_cache', flag)

    def set_cache_max_size(self, diff):
        parms.set_int_par(
                self.handle, self.name, 'cache_max_size', diff)

    def get_cache_max_size(self):
        return parms.int_par(self.handle, self.name, 'cache_max_size')

    def set_subsensitivity_filenames(self, names):
        return parms.set_char_par(self.handle, self.name, 'subsensitivity_filenames', names)

    def get_subsensitivity_filenames(self):
        return parms.char_par(self.handle, self.name, 'subsensitivity_filenames')


class Reconstructor(object):
    """Base class for a generic PET reconstructor."""

    def __init__(self):
        """init."""
        self.handle = None
        self.input = None
        self.image = None

    def __del__(self):
        """del."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_input(self, input_data):
        """Sets the acquisition data to use for reconstruction."""
        assert_validity(input_data, ScanData)
        parms.set_parameter(
            self.handle, 'Reconstruction',
            'input_data', input_data.handle)

    def set_output_filename_prefix(self, prefix):
        """Sets the output file name prefix."""
        parms.set_char_par(
            self.handle, 'Reconstruction', 'output_filename_prefix', prefix)

    def disable_output(self):
        """Disables output."""
        parms.set_int_par(self.handle, 'Reconstruction', 'disable_output', 1)

    def enable_output(self):
        """Enables output."""
        parms.set_int_par(self.handle, 'Reconstruction', 'enable_output', 1)

    def reconstruct(self, image):
        """Performs reconstruction (will update the image argument)"""
        assert_validity(image, ImageData)
        try_calling(pystir.cSTIR_runReconstruction(self.handle, image.handle))
        self.image = image

    def get_output(self):
        """Returns the reconstructed image."""
        assert self.image is not None, 'current estimate not set. Did you run a reconstruction already?'
        # TODO: move to C++
        return self.image.clone()


class FBP2DReconstructor(object):
    """Class for 2D Filtered Back Projection reconstructor.

    This is an implementation of the 2D FBP algorithm.
    Oblique angles in data will be ignored. The exception is the span=1 case,
    where the ring differences +1 and -1 are first combined to give indirect
    sinograms.
    By default, the algorithm uses the ramp filter. An apodizing filter can be
    added by using set_alpha_cosine_window and/or set_frequency_cut_off.
    The apodizing filter in frequency space has the form

        (alpha + (1 - alpha) * cos(pi * f / fc))
    """

    def __init__(self):
        """init."""
        self.handle = None
        self.handle = pystir.cSTIR_newObject('FBP2D')
        check_status(self.handle)

    def __del__(self):
        """del."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_input(self, input_data):
        """Sets the acquisition data to use for reconstruction."""
        assert_validity(input_data, AcquisitionData)
        parms.set_parameter(self.handle, 'FBP2D', 'input', input_data.handle)

    def set_zoom(self, v):
        """Sets zoom."""
        parms.set_float_par(self.handle, 'FBP2D', 'zoom', v)

    def set_alpha_cosine_window(self, v):
        """Sets alpha in the apodizing filter.

        See the class documentation for the filter. The value of alpha should
        be between 0.5 and 1. alpha=0.5 corresponds to the Hann filter, while
        0.54 corresponds to the Hamming filter.
        """
        parms.set_float_par(self.handle, 'FBP2D', 'alpha', v)

    def set_frequency_cut_off(self, v):
        """Sets the cut-off frequency for the apodizing filter.

        See the class documentation for the filter. The value of fc should be
        between 0 and 0.5.
        """
        parms.set_float_par(self.handle, 'FBP2D', 'fc', v)

    def set_output_image_size_xy(self, xy):
        """Sets output image size (xy)."""
        parms.set_int_par(self.handle, 'FBP2D', 'xy', xy)

    def set_up(self, image):
        """Sets up the reconstructor."""
        try_calling(pystir.cSTIR_setupFBP2DReconstruction(
            self.handle, image.handle))

    def process(self):
        """Performs reconstruction."""
        try_calling(pystir.cSTIR_runFBP2DReconstruction(self.handle))

    def get_output(self):
        """Returns the reconstructed image."""
        image = ImageData()
        image.handle = parms.parameter_handle(self.handle, 'FBP2D', 'output')
        check_status(image.handle)
        return image


class IterativeReconstructor(Reconstructor):
    """Base class for a generic iterative PET reconstructor."""

    def __init__(self):
        """init."""
        self.handle = None
        self.image = None
        self.subset = 0

    def __del__(self):
        """del."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_num_subsets(self, n):
        """See method in ObjectiveFunction."""
        parms.set_int_par(
            self.handle, 'IterativeReconstruction', 'num_subsets', n)

    def get_num_subsets(self):
        """Returns number of subsets."""
        return parms.int_par(
            self.handle, 'IterativeReconstruction', 'num_subsets')
#    def set_start_subset_num(self, n):
#        parms.set_int_par\
#            (self.handle, 'IterativeReconstruction', 'start_subset_num', n)
#    def get_start_subset_num(self):
#        return parms.int_par\
#            (self.handle, 'IterativeReconstruction', 'start_subset_num')

    def set_num_subiterations(self, n):
        """Sets number of subiterations.

        In subset iterative methods, each iteration works with a subset,
        and is therefore referred to as a subiteration.
        """
        parms.set_int_par(
            self.handle, 'IterativeReconstruction', 'num_subiterations', n)

    def get_num_subiterations(self):
        """Returns the number of subiterations."""
        return parms.int_par(
            self.handle, 'IterativeReconstruction', 'num_subiterations')
#    def set_start_subiteration_num(self, n):
#        parms.set_int_par(
#           self.handle, 'IterativeReconstruction',
#           'start_subiteration_num', n)
#    def get_start_subiteration_num(self):
#        return parms.int_par\
#            (self.handle, 'IterativeReconstruction', 'start_subiteration_num')
#    def set_subiteration_num(self, iter):
#        parms.set_int_par\
#            (self.handle, 'IterativeReconstruction', 'subiteration_num', iter)

    def get_subiteration_num(self):
        """Returns the current subiteration number."""
        return parms.int_par(
            self.handle, 'IterativeReconstruction', 'subiteration_num')

    def set_save_interval(self, n):
        """Sets save interval.

        Defines how often to save image iterates (n = 1: on each
        subiteration, n = 2: every other subiteration etc.)
        """
        if n > 0:
            self.enable_output()
        else:
            self.disable_output()
        parms.set_int_par(
            self.handle, 'IterativeReconstruction', 'save_interval', n)
#    def set_inter_iteration_filter_interval(self, n):
#        parms.set_int_par\
#            (self.handle, 'IterativeReconstruction',\
#             'inter_iteration_filter_interval', n)

    def set_objective_function(self, obj):
        """Specifies the objective function to be maximized."""
        assert_validity(obj, ObjectiveFunction)
        parms.set_parameter(self.handle, 'IterativeReconstruction',
            'objective_function', obj.handle)

    abc.abstractmethod
    def get_objective_function(self):
        pass
#        obj_fun = ObjectiveFunction()
#        obj_fun.handle = pystir.cSTIR_parameter\
#            (self.handle, 'IterativeReconstruction', 'objective_function')
#        check_status(obj_fun.handle)
#        return obj_fun
#    def set_inter_iteration_filter(self, f):
#        pystir.cSTIR_setParameter\
#            (self.handle, 'IterativeReconstruction',\
#             'inter_iteration_filter_type', f.handle)
#    def get_inter_iteration_filter(self):
#        filter = DataProcessor()
#        filter.handle = pystir.cSTIR_parameter\
#            (self.handle, 'IterativeReconstruction',\
#             'inter_iteration_filter_type')
#        check_status(filter.handle)
#        return filter

    def set_up(self, image):
        """Sets up the reconstructor."""
        assert_validity(image, ImageData)
        try_calling(pystir.cSTIR_setupReconstruction(
            self.handle, image.handle))

    def set_current_estimate(self, image):
        """Sets image estimate for further iterations.

        image will be cloned.
        """
        assert_validity(image, ImageData)
        self.image = image.clone()

    #def set_estimate(self, image):
    #    """Sets image estimate as a variable that will be updated."""
    #    assert_validity(image, ImageData)
    #    self.image = image

    def process(self):
        """Performs reconstruction."""
        assert self.image is not None, 'current estimate not set.'
        try_calling(pystir.cSTIR_runReconstruction(
            self.handle, self.image.handle))

    def get_current_estimate(self):
        """Return current image estimate."""
        assert self.image is not None, 'current estimate not set.'
        return self.image.clone()

    def update_current_estimate(self):
        """Updates current image estimate by performing one subiteration."""
        if self.image is None:
            raise error('current estimate not set')
        assert_validity(self.image, ImageData)
        try_calling(pystir.cSTIR_updateReconstruction(
            self.handle, self.image.handle))

    def set_current_subset_num(self, subset):
        """Sets the subset to be used on the next subiteration."""
        self.subset = subset

    def get_subset_sensitivity(self):
        """Returns the current subset sensitivity.

        Returns an ImageData object containing sensitivity image for the
        current subset.
        """
        obj_fun = self.get_objective_function()
        return obj_fun.get_subset_sensitivity(self.subset)

    def update(self, image):
        """Updates the image estimate.

        Applies one subiteration to the image estimate passed as the
        argument.
        """
        assert_validity(image, ImageData)
        #self.set_estimate(image)
        self.image = image;
        self.update_current_estimate()
        return image


class OSMAPOSLReconstructor(IterativeReconstructor):
    """OSMAPOSL reconstruction class.

    Class for reconstructor objects using Ordered Subsets Maximum A Posteriori
    One Step Late reconstruction algorithm, see
    http://stir.sourceforge.net/documentation/doxy/html/classstir_1_1OSMAPOSLReconstruction.html
    """

    def __init__(self, filename=''):
        """init."""
        self.handle = None
        self.image = None
        self.name = 'OSMAPOSL'
        self.handle = pystir.cSTIR_objectFromFile(
            'OSMAPOSLReconstruction', filename)
        check_status(self.handle)
        self.disable_output()

    def __del__(self):
        """del."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_maximum_relative_change(self, value):
        """Sets maximum relative change."""
        parms.set_float_par(
            self.handle, self.name, 'set_maximum_relative_change', value)

    def set_minimum_relative_change(self, value):
        """Sets minimum relative change."""
        parms.set_float_par(
            self.handle, self.name, 'set_minimum_relative_change', value)
#    def set_MAP_model(self, model):
#        parms.set_char_par\
#            (self.handle, self.name, 'MAP_model', model)
    def get_objective_function(self):
        obj_fun = PoissonLogLikelihoodWithLinearModelForMean()
        obj_fun.handle = pystir.cSTIR_parameter\
            (self.handle, self.name, 'objective_function')
        check_status(obj_fun.handle)
        return obj_fun


class KOSMAPOSLReconstructor(IterativeReconstructor):
    """KOSMAPOSL recontstructor class.

    Class for reconstructor objects using Kernel Ordered Subsets Maximum
    A Posteriori One Step Late reconstruction algorithm.

    This class implements the iterative algorithm obtained using the Kernel
    method (KEM) and Hybrid kernel method (HKEM). This implementation
    corresponds to the one presented by Deidda D et al, "Hybrid PET-MR
    list-mode kernelized expectation maximization  reconstruction", Inverse
    Problems, 2019, DOI: https://doi.org/10.1088/1361-6420/ab013f.
    However, this allows also sinogram-based reconstruction. Each voxel value
    of the image X can be represented as a linear combination using the kernel
    method.  If we have an image with prior information, we can construct for
    each voxel j of the emission image a feature vector, v, using the prior
    information. The image X can then be described using the kernel matrix

    X = A*K

    where K is the kernel matrix. The resulting algorithm with OSEM,
    for example, is the following:

    A^(n+1) =  A^n/(K^n * S) * K^n * P * Y/(P * K^n *A^n + S)

    where kernel can be written as:

    K^n = K_m * K_p;

    with

    K_m = exp(-(v_j - v_l)^2/(2*sigma_m^2)) *
          exp(-(x_j - x_l)^2 /(2*sigma_dm^2))

    being the MR component of the kernel and

    K_p = exp(-(z_j - z_l)^2/(2*sigma_p^2)) *
          exp(-(x_j - x_l)^2 /(2*sigma_dp^2))

    is the part coming from the emission iterative update. Here, the Gaussian
    kernel functions have been modulated by the distance between voxels in the
    image space.
    """

    def __init__(self, filename=''):
        """init."""
        IterativeReconstructor.__init__(self)
        self.handle = None
        self.image = None
        self.name = 'KOSMAPOSL'
        self.handle = pystir.cSTIR_objectFromFile(
            'KOSMAPOSLReconstruction', filename)
        check_status(self.handle)
        self.disable_output()

    def __del__(self):
        """del."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_anatomical_prior(self, ap):
        """Sets anatomical prior."""
        assert_validity(ap, ImageData)
        parms.set_parameter(
            self.handle, 'KOSMAPOSL', 'anatomical_prior', ap.handle)

    def set_num_neighbours(self, n):
        """Sets number of neighbours."""
        parms.set_int_par(
            self.handle, 'KOSMAPOSL', 'num_neighbours', n)

    def set_num_non_zero_features(self, n):
        """Sets number of non-zero features."""
        parms.set_int_par(
            self.handle, 'KOSMAPOSL', 'num_non_zero_features', n)

    def set_sigma_m(self, v):
        """Sets sigma m."""
        parms.set_float_par(self.handle, 'KOSMAPOSL', 'sigma_m', v)

    def set_sigma_p(self, v):
        """Sets sigma p."""
        parms.set_float_par(self.handle, 'KOSMAPOSL', 'sigma_p', v)

    def set_sigma_dm(self, v):
        """Sets sigma dm."""
        parms.set_float_par(self.handle, 'KOSMAPOSL', 'sigma_dm', v)

    def set_sigma_dp(self, v):
        """Sets sigma dp."""
        parms.set_float_par(self.handle, 'KOSMAPOSL', 'sigma_dp', v)

    def set_only_2D(self, tf):
        """Sets only_2D flag."""
        v = 1 if tf else 0
        parms.set_int_par(self.handle, 'KOSMAPOSL', 'only_2D', v)

    def set_hybrid(self, tf):
        """Sets use hybrid mode flag."""
        v = 1 if tf else 0
        parms.set_int_par(self.handle, 'KOSMAPOSL', 'hybrid', v)

    def compute_kernelised_image(self, image, alpha):
        assert_validity(image, ImageData)
        assert_validity(alpha, ImageData)
        ki = ImageData()
        ki.handle = pystir.cSTIR_computeKernelisedImage \
            (self.handle, image.handle, alpha.handle)
        check_status(ki.handle)
        return ki

    def get_objective_function(self):
        obj_fun = PoissonLogLikelihoodWithLinearModelForMean()
        obj_fun.handle = pystir.cSTIR_parameter\
            (self.handle, self.name, 'objective_function')
        check_status(obj_fun.handle)
        return obj_fun

class SingleScatterSimulator():
    '''
    Class for simulating the scatter contribution to PET data.

    This class uses the STIR Single Scatter simulation, taking as input an
    activity and attenuation image, and a acquisition data template.

    WARNING: Currently this class does not use the low-resolution sampling
    mechanism of STIR. This means that if you give it a full resolution acq_data,
    you will likely run out of memory and/or time.
    '''
    def __init__(self, filename = ''):
        self.handle = None
        self.image = None
        self.name = 'PETSingleScatterSimulator'
        self.filename = filename

        if not self.filename:
            self.handle = pystir.cSTIR_newObject(self.name)
        else:
            self.handle = pystir.cSTIR_objectFromFile(self.name, self.filename)
        check_status(self.handle)

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_up(self, acq_templ, img_templ):
        """Set up.

        Prepare this object for performing forward operations;
        acq_templ:  an AcquisitionData object used as a template for
                    creating an AcquisitionData object to store forward
                    projection;
        img_templ:  an ImageData object used as a template for checking geometry etc

        attenuation image has to be set first
        """
        assert_validity(acq_templ, AcquisitionData)
        assert_validity(img_templ, ImageData)

        # temporarily save the templates in the class
        self.acq_templ = acq_templ
        #self.img_templ = img_templ

        try_calling(pystir.cSTIR_setupScatterSimulator(
            self.handle, acq_templ.handle, img_templ.handle))

    def forward(self, image,  out=None):
        """Return the scatter estimation for the input activity image.

        image   :  an ImageData object.

        set_up() has to be called first.
        """
        assert_validity(image, ImageData)
        if out is None:
            ad = AcquisitionData()
            ad.handle = pystir.cSTIR_scatterSimulatorFwd(
                self.handle, image.handle);
            check_status(ad.handle)
            return ad
        ad = out
        assert_validity(ad, AcquisitionData)
        try_calling(pystir.cSTIR_scatterSimulatorFwdReplace(
            self.handle, image.handle, ad.handle))

    def set_attenuation_image(self, image):
        assert_validity(image, ImageData)
        parms.set_parameter(self.handle, self.name, 'setAttenuationImage', image.handle)


class ScatterEstimator():
    '''
    Class for estimating the scatter contribution in PET projection data

    This class implements the SSS iterative algorithm from STIR. It
    is an iterative loop of reconstruction, single scatter estimation,
    upsampling, tail-fitting.

    Output is an acquisition_data object with the scatter contribution.
    This can then be added to the randoms to use in PETAcquisitionModel.set_background_term().
    '''
    def __init__(self, filename = ''):
        self.handle = None
        self.image = None
        self.name = 'PETScatterEstimator'
        self.filename = filename

        if not self.filename:
            self.handle = pystir.cSTIR_newObject(self.name)
        else:
            self.handle = pystir.cSTIR_objectFromFile(self.name, self.filename)

        check_status(self.handle)

    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_up(self):
        """
        Set up.

        Prepare this object for performing scatter estimation;
        All input has to be set before calling this function.
        """
        try_calling(pystir.cSTIR_setupScatterEstimator(
            self.handle))

    def process(self):
        """
        Runs the scatter estimation.

        You need to run set_up() first.
        """
        print('ScatterEstimator:: Waiting for the scatter estimation to finish ...')
        self.output = AcquisitionData()
        self.output.handle = pystir.cSTIR_runScatterEstimator(self.handle)
        check_status(self.output.handle)
        print('ScatterEstimator:: estimation finished.')

    def get_output(self):
        """
        Return the final scatter estimate.
        """
        data = AcquisitionData()
        data.handle = parms.parameter_handle(self.handle, 'PETScatterEstimator', 'output')
        check_status(data.handle)
        return data

    def get_num_iterations(self):
        """Get number of iterations of the SSS algorithm to use."""
        return parms.int_par(self.handle, 'PETScatterEstimator', 'num_iterations')

    def get_OSEM_num_subiterations(self):
        """Get number of subiterations used by OSEM in the SSS algorithm."""
        return parms.int_par(self.handle, 'PETScatterEstimator', 'OSEM_num_subiterations')

    def get_OSEM_num_subsets(self):
        """Get number of subsets used by OSEM in the SSS algorithm."""
        return parms.int_par(self.handle, 'PETScatterEstimator', 'OSEM_num_subsets')

    def set_attenuation_image(self, image):
        assert_validity(image, ImageData)
        parms.set_parameter(self.handle, self.name, 'setAttenuationImage', image.handle)

    def set_attenuation_correction_factors(self, arg):
        assert_validity(arg, AcquisitionData)
        parms.set_parameter(self.handle, self.name, 'setAttenuationCorrectionFactors', arg.handle)

    def set_mask_image(self, image):
        assert_validity(image, ImageData)
        parms.set_parameter(self.handle, self.name, 'setMaskImage', image.handle)

    def set_mask_acq_data(self, arg):
        assert_validity(arg, AcquisitionData)
        parms.set_parameter(self.handle, self.name, 'setMaskAcqData', arg.handle)


    def set_input(self, acq_data):
        assert_validity(acq_data, AcquisitionData)
        parms.set_parameter(self.handle, self.name, 'setInput', acq_data.handle)

    def set_randoms(self, acq_data):
        assert_validity(acq_data, AcquisitionData)
        parms.set_parameter(self.handle, self.name, 'setRandoms', acq_data.handle)

    def set_asm(self, asm):
        '''Set acquisition sensitivity model (without attenuation!)'''
        assert_validity(asm, AcquisitionSensitivityModel)
        parms.set_parameter(self.handle, self.name, 'setASM', asm.handle)

    def set_OSEM_num_subiterations(self, v):
        """Set number of subiterations used by OSEM in the SSS algorithm."""
        parms.set_int_par(self.handle, 'PETScatterEstimator', 'set_OSEM_num_subiterations', v)

    def set_OSEM_num_subsets(self, v):
        """Set number of subsets used by OSEM in the SSS algorithm."""
        parms.set_int_par(self.handle, 'PETScatterEstimator', 'set_OSEM_num_subsets', v)

    def set_num_iterations(self, v):
        """Set number of iterations of the SSS algorithm to use."""
        parms.set_int_par(self.handle, 'PETScatterEstimator', 'set_num_iterations', v)

    def set_max_scale_value(self, v):
        """Set maximal scale factor value of the SSS algorithm to use."""
        parms.set_float_par(self.handle, 'PETScatterEstimator', 'set_max_scale_value', v)

    def set_min_scale_value(self, v):
        """Set maximal scale factor value of the SSS algorithm to use."""
        parms.set_float_par(self.handle, 'PETScatterEstimator', 'set_min_scale_value', v)

    def set_output_prefix(self, v):
        """
        Set prefix for filenames with scatter estimates.

        Actual filenames will append the iteration number and the .hs extension
        as common for STIR Interfile data.

        Set it to the empty string to prevent any output.
        """
        parms.set_char_par(self.handle, 'PETScatterEstimator', 'set_output_prefix', v)

class OSSPSReconstructor(IterativeReconstructor):
    """OSSPS reconstructor class.

    Class for reconstructor objects using Ordered Subsets Separable
    Paraboloidal Surrogate reconstruction algorithm, see
    http://stir.sourceforge.net/documentation/doxy/html/classstir_1_1OSSPSReconstruction.html
    """

    def __init__(self, filename=''):
        """init."""
        self.handle = None
        self.image = None
        self.name = 'OSSPS'
        self.handle = pystir.cSTIR_objectFromFile(
            'OSSPSReconstruction', filename)
        check_status(self.handle)
        self.disable_output()

    def __del__(self):
        """del."""
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

    def set_relaxation_parameter(self, value):
        """Sets relaxation parameter."""
        parms.set_float_par(
            self.handle, self.name, 'relaxation_parameter', value)

    def set_relaxation_gamma(self, value):
        """Sets relaxation gamma parameter."""
        parms.set_float_par(
            self.handle, self.name, 'relaxation_gamma', value)

    def set_upper_bound(self, value):
        """Sets upper bound parameter."""
        parms.set_double_par(
            self.handle, self.name, 'upper_bound', value)

    def get_relaxation_parameter(self):
        """Returns relaxation parameter value."""
        return parms.float_par(self.handle, self.name, 'relaxation_parameter')

    def get_relaxation_gamma(self):
        """Returns relaxation gamma value."""
        return parms.float_par(self.handle, self.name, 'relaxation_gamma')

    def get_upper_bound(self):
        """Returns upper bound value."""
        return parms.double_par(self.handle, self.name, 'upper_bound')


#def make_Poisson_loglikelihood(acq_data, likelihood_type='LinearModelForMean',
def make_Poisson_loglikelihood(acq_data=None, likelihood_type=None,
                               acq_model=None):
    """Makes Poisson loglikelihood.

    Selects the objective function based on the acquisition data and
    one of the following likelihood model types:
    * LinearModelForMean
    * LinearModelForMeanAndListModeDataWithProjMatrixByBin
    """
    if likelihood_type is None or likelihood_type=='LinearModelForMean':
        obj_fun = PoissonLogLikelihoodWithLinearModelForMeanAndProjData()
        if acq_data is not None:
            obj_fun.set_acquisition_data(acq_data)
        else:
            raise error('PoissonLogLikelihoodWithLinearModelForMeanAndProjData' + \
            ' requires acquisition data')
    elif likelihood_type == 'LinearModelForMeanAndListModeDataWithProjMatrixByBin':
        obj_fun = PoissonLogLikelihoodWithLinearModelForMeanAndListModeDataWithProjMatrixByBin()
        if acq_data is not None:
            obj_fun.set_acquisition_data(acq_data)
    else:
        raise error('Poisson_loglikelihood of type ' + likelihood_type + \
                    ' is not implemented')
    if acq_model is not None:
        obj_fun.set_acquisition_model(acq_model)
    return obj_fun
