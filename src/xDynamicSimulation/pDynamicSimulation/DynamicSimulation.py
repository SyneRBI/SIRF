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
import numpy as np
from numbers import Number, Complex, Integral


try:
    import pylab
    HAVE_PYLAB = True
except:
    HAVE_PYLAB = False

import sys

if sys.version_info[0] >= 3 and sys.version_info[1] >= 4:
    ABC = abc.ABC
else:
    ABC = abc.ABCMeta('ABC', (), {})

from deprecation import deprecated

from sirf.Utilities import check_status, assert_validity, try_calling, error

import sirf
from sirf import SIRF
from sirf.SIRF import DataContainer

import sirf.pysimulation as pysim
import sirf.pyiutilities as pyiutil
# import sirf.pygadgetron as pygadgetron
# import sirf.pysirf as pysirf

import sirf.Reg as pReg
import sirf.Gadgetron as pMR



def format_arrays_for_setters(data):

    if not isinstance(data, np.ndarray):
        raise error('Wrong input format.' + \
            ' Should be numpy.ndarray. Got {}'.format(type(data)))

    if data.dtype != np.float32:
            the_data = data.astype(np.float32)

    convert = not the_data.flags['C_CONTIGUOUS']
    if convert:
        the_data = np.ascontiguousarray(the_data)

    return the_data


class DynamicSimulation(object):

    def __init__(self, tissue_labels, fname_xml):

        self.handle = None
        assert_validity(tissue_labels, pReg.NiftiImageData3D)

        self.handle = pysim.cDS_MRDynamicSimulation(tissue_labels.handle, fname_xml)
        check_status(self.handle)

    def set_acquisition_template_data(self, ad):
        assert_validity(ad, pMR.AcquisitionData)
        pysim.cDS_setAcquisitionTemplateData(self.handle, ad.handle)

    def simulate_data(self):
        pysim.cDS_simulateData(self.handle)

    def set_csm(self, csm):
        assert_validity(csm, pMR.CoilSensitivityData)
        pysim.cDS_setCoilmaps(self.handle, csm.handle)

    def set_snr(self, SNR):
        pysim.cDS_setSNR(self.handle, SNR)

    def set_snr_label(self, SNR_label):
        pysim.cDS_setNoiseLabel(self.handle, SNR_label)

    def write_simulation_results(self, fpath_output):
        pysim.cDS_writeSimulationResults(self.handle, fpath_output); 

    def add_motion_dynamic(self, motiondyn):
        pysim.cDS_addMRMotionDynamic(self.handle, motiondyn.handle); 



class SurrogateSignal(object):
    def __init__(self, time_points, signal_points):

        if not isinstance(time_points, np.ndarray):
            raise error('Wrong input format.' + \
                ' Should be numpy.ndarray. Got {}'.format(type(time_points)))

        if not isinstance(signal_points, np.ndarray):
            raise error('Wrong input format.' + \
                ' Should be numpy.ndarray. Got {}'.format(type(signal_points)))

        num_signal_points = time_points.size
        if num_signal_points != signal_points.size:
            AssertionError("The signal and time points do not have the identical size.")

        if time_points.dtype != np.float32:
            the_time = time_points.astype(np.float32)
        
        if signal_points.dtype != np.float32:
            the_signal = signal_points.astype(np.float32)

        convert = not the_time.flags['C_CONTIGUOUS']
        if convert:
            the_time = np.ascontiguousarray(the_time)

        convert = not the_signal.flags['C_CONTIGUOUS']
        if convert:
            the_signal = np.ascontiguousarray(the_signal)
        
        self.handle = None
        self.handle = pysim.cDS_DynamicSignal(the_time.ctypes.data,\
                                    the_signal.ctypes.data,\
                                    num_signal_points)
        check_status(self.handle)


class Dynamic(object):

    def set_dynamic_signal(self, time_points, signal_points):
      
        num_signal_points = time_points.size
        if num_signal_points != signal_points.size:
            AssertionError("The signal and time points do not have the identical size.")

        the_time = format_arrays_for_setters(time_points)
        the_signal = format_arrays_for_setters(signal_points)
        
        pysim.cDS_setDynamicSignal(self.handle, 
                                    the_time.ctypes.data,\
                                    the_signal.ctypes.data,\
                                    num_signal_points)

    def add_displacement_field(self, dvf):
        assert_validity(dvf, pReg.NiftiImageData3DDisplacement)
        pysim.cDS_addMRDisplacementField(self.handle, dvf.handle)

    def set_cyclicality(self, is_cyclic):
        pysim.cDS_setCyclicality(self.handle, is_cyclic)

    def get_num_signal_pts(self):
        
        h_pts = pysim.cDS_getNumberOfSignalPoints(self.handle)
        return pyiutil.intDataFromHandle(h_pts)

class MRDynamic(Dynamic):
    
    def set_mr_acquisitions(self, ad):
        assert_validity(ad, pMR.AcquisitionData)
        pysim.cDS_setMRAcquisitions(self.handle, ad.handle)


class MRMotionDynamic(MRDynamic):

    def __init__(self, num_states):
        self.handle = None
        self.handle = pysim.cDS_MRMotionDynamic(num_states)
        check_status(self.handle)
   
    def add_displacement_field(self, dvf):
        assert_validity(dvf, pReg.NiftiImageData3DDisplacement)
        pysim.cDS_addMRDisplacementField(self.handle, dvf.handle)

    def set_cyclicality(self, is_cyclic):
        pysim.cDS_setCyclicality(self.handle, is_cyclic)

# #class ImageData(DataContainer):
# class ImageData(SIRF.ImageData):
#     '''
#     Class for an MR images container.
#     Each item in the container is a 3D complex or float array of the image 
#     values on an xyz-slice (z-dimension is normally 1).
#     '''
#     def __init__(self, file = None):
#         self.handle = None
#         if file is None:
#             return
#         self.handle = pygadgetron.cGT_readImages(file)
#         check_status(self.handle)
#     def __del__(self):
#         if self.handle is not None:
#             pyiutil.deleteDataHandle(self.handle)
#     def same_object(self):
#         return ImageData()

#     def abs(self):
#         print('in ImageData.abs()...')
#         images = ImageData()
#         images.handle = pygadgetron.cGT_absImages(self.handle)
#         check_status(images.handle)
#         return images

#     def read_from_file(self, file):
#         if self.handle is not None:
#             pyiutil.deleteDataHandle(self.handle)
#         self.handle = pygadgetron.cGT_readImages(file)
#         check_status(self.handle)

#     def from_acquisition_data(self, ad):
#         assert isinstance(ad, AcquisitionData), "Please pass a AcquisitionData object"
#         if self.handle is not None:
#             pyiutil.deleteDataHandle(self.handle)
#         self.handle = pygadgetron.cGT_ImageFromAcquisitiondata(ad.handle)
#         check_status(self.handle)

#     def data_type(self, im_num):
#         '''
#         Returns the data type for a specified image (see 8 data types above).
#         im_num: image (slice) 
#         '''
#         assert self.handle is not None
#         handle = pygadgetron.cGT_imageDataType(self.handle, im_num)
#         check_status(handle)
#         n = pyiutil.intDataFromHandle(handle)
#         pyiutil.deleteDataHandle(handle)
#         return n
#     def is_real(self):
#         assert self.handle is not None
#         t = self.data_type(0)
#         return t is not ISMRMRD_CXFLOAT and t is not ISMRMRD_CXDOUBLE
#     def process(self, list):
#         '''
#         Returns processed self with an image processor specified by
#         a list of gadgets.
#         list: Python list of gadget description strings, each gadget 
#               description being a string of the form
#                 '[label:]gadget_name[(property1=value1[,...])]'
#               (square brackets embrace optional items, ... stands for etc.)
#         '''
#         assert self.handle is not None
#         ip = ImageDataProcessor(list)
#         return ip.process(self)
#     def image(self, im_num):
#         return Image(self, im_num)
#     def select(self, attr, value):
#         '''
#         Creates an images container with images from self with the specified
#         value of specified attribute.
#         attr : the name of the attribute (Python string)
#         value: the value of the attribute (Python string)
#         '''
#         assert self.handle is not None
#         images = ImageData()
#         images.handle = pygadgetron.cGT_selectImages(self.handle, attr, value)
#         check_status(images.handle)
#         return images

#     def get_ISMRMRD_info(self, par):
#         '''
#         Returns the array of values of the specified image information 
#         parameter. Parameters names are the same as the names of sirf.Gadgetron.Image class
#         public methods (except is_real and info).

#         par: parameter name (as a string)

#         Examples:

#         # to get information on the timing:
#         slice = image.get_ISMRMRD_info('acquisition_time_stamp')

#         # to get the unit vector orthogonal to the slice and directed
#         # to the next slice (in ISMRMRD coordinate system):
#         slice_dir = image.get_ISMRMRD_info('slice_dir')
#         '''
#         ni = self.number()
#         info = numpy.empty((ni,), dtype = object)
#         for i in range(ni):
#             image = self.image(i)
#             info[i] = image.info(par)
#         return info

#     @deprecated(details="Please use get_ISMRMRD_info method instead")
#     def get_info(self, par):
#         return self.get_ISMRMRD_info(par)

#     def fill(self, data):
#         '''
#         Fills self's image data with specified values.
#         data: Python Numpy array or ImageData
#         '''
#         assert self.handle is not None
        
#         if isinstance(data, ImageData):
#             return super(ImageData, self).fill(data)
        
#         if isinstance(data, numpy.ndarray):
#             the_data = data
#             if self.is_real():
#                 if data.dtype != numpy.float32:
#                     the_data = data.astype(numpy.float32)
#             else:
#                 if data.dtype != numpy.complex64:
#                     the_data = data.astype(numpy.complex64)
#             convert = not data.flags['C_CONTIGUOUS']
#             if convert:
#                 the_data = numpy.ascontiguousarray(the_data)
#             if self.is_real():
#                 try_calling(pygadgetron.cGT_setImageDataFromFloatArray\
#                     (self.handle, the_data.ctypes.data))
#             else:
#                 try_calling(pygadgetron.cGT_setImageDataFromCmplxArray\
#                     (self.handle, the_data.ctypes.data))
#         elif isinstance (data, Complex):
#             arr = data + numpy.zeros(self.shape, dtype=numpy.complex64)
#             return self.fill(arr)
#         elif isinstance (data, Number):
#             arr = data + numpy.zeros(self.shape, dtype=numpy.float32)
#             return self.fill(arr)
#         else:
#             raise error('wrong fill value.' + \
#                         ' Should be ImageData, numpy.ndarray or number. Got {}'.format(type(data)))
#         return self

#     def dimensions(self):
#         '''
#         Returns the dimensions of 3D/4D Numpy ndarray of all self's images.
#         '''
#         if self.number() < 1:
#             return 0
#         assert self.handle is not None
#         dim = numpy.ndarray((4,), dtype = numpy.int32)
#         image = Image(self)
#         pygadgetron.cGT_getImageDim(image.handle, dim.ctypes.data)
#         nx = dim[0]
#         ny = dim[1]
#         nz = dim[2]
#         nc = dim[3]
        
#         nz = nz*self.number()
        
#         if nc == 1: # for backward compatibility
#             return nz, ny, nx
#         else:
#             return nc, nz, ny, nx
        
#     def as_array(self):
#         '''
#         Returns all self's images as a 3D or 4D Numpy ndarray.
#         '''
#         dims = self.dimensions()
        
#         assert self.handle is not None
#         if self.number() < 1:
#             return numpy.ndarray((0,), dtype = numpy.float32)
#         if self.is_real():
#             array = numpy.ndarray(dims, dtype = numpy.float32)
#             try_calling(pygadgetron.cGT_getImageDataAsFloatArray\
#                 (self.handle, array.ctypes.data))
#         else:
#             array = numpy.ndarray(dims, dtype = numpy.complex64)
#             try_calling(pygadgetron.cGT_getImageDataAsCmplxArray\
#                 (self.handle, array.ctypes.data))
                
#         if len(dims) != 4:
#             return array

#         nc, nz, ny, nx = dims
#         ns = self.number() # number of total dynamics (slices, contrasts, etc.)
#         nz = nz//ns        # z-dimension of a slice

#         # hope Numpy is clever enough to do all this in-place:
#         array = numpy.reshape(array, (ns, nc, nz, ny, nx))
#         array = numpy.swapaxes(array, 0, 1)
#         array = numpy.reshape(array, (nc, ns*nz, ny, nx))

#         return array
                
#     def copy(self):
#         '''alias of clone'''
#         return self.clone()
#     def conjugate(self):
#         '''Returns the complex conjugate of the data '''
#         if self.handle is not None:
#             out = self.clone()
#             out.fill(self.as_array().conjugate())
#             return out
#         else:
#             raise error("Empty object cannot be conjugated")
#     def show(self, zyx=None, slice=None, title=None, cmap='gray', postpone=False):
#         '''Displays xy-cross-section(s) of images.'''
#         assert self.handle is not None
#         if not HAVE_PYLAB:
#             print('pylab not found')
#             return
#         data = self.as_array()
#         nz = data.shape[0]
#         if isinstance(slice, (Integral,numpy.integer)):
#             if slice < 0 or slice >= nz:
#                 return
#             ni = 1
#             slice = [slice]
#         elif slice is None:
#             ni = nz
#             slice = range(nz)
#         else:
#             try:
#                 ni = len(slice)
#             except:
#                 raise error('wrong slice list')
#         if title is None:
#             title = 'Selected images'
#         if ni >= 16:
#             tiles = (4, 4)
#         else:
#             tiles = None
#         f = 0
#         while f < ni:
#             t = min(f + 16, ni)
#             err = show_3D_array(abs(data), index=slice[f : t], \
#                                 tile_shape=tiles, cmap=cmap, \
#                                 zyx=zyx, label='image', \
#                                 xlabel='samples', ylabel='readouts', \
#                                 suptitle=title, \
#                                 show=(t == ni) and not postpone)
#             f = t
#     def allocate(self, value=0, **kwargs):
#         '''Method to allocate an ImageData and set its values
        
#         CIL/SIRF compatibility
#         '''
#         if value in ['random', 'random_int']:
#             out = self.clone()
#             shape = out.as_array().shape
#             seed = kwargs.get('seed', None)
#             if seed is not None:
#                 numpy.random.seed(seed) 
#             if value == 'random':
#                 out.fill(numpy.random.random_sample(shape))
#             elif value == 'random_int':
#                 max_value = kwargs.get('max_value', 100)
#                 out.fill(numpy.random.randint(max_value,size=shape))
#         else:
#             out = self.clone()
#             tmp = value * numpy.ones(out.as_array().shape)
#             out.fill(tmp)
#         return out

#     def print_header(self, im_num):
#         """Print the header of one of the images. zero based."""
#         try_calling(pygadgetron.cGT_print_header(self.handle, im_num))
#     @property
#     def dtype(self):
#         if self.is_real():
#             return numpy.float32
#         return numpy.complex64
#     @property
#     def shape(self):
#         return self.dimensions()
        
# SIRF.ImageData.register(ImageData)



