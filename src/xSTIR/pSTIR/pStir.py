''' Object-Oriented wrap for low-level STIR Python interface pystir.py
'''
import numpy
import os
try:
    import pylab
    HAVE_PYLAB = True
except:
    HAVE_PYLAB = False
import sys
import time

from pUtil import *
import pyiutil
import pystir

INFO_CHANNEL = 0
WARNING_CHANNEL = 1
ERROR_CHANNEL = 2
ALL_CHANNELS = -1

def pet_data_path():
    return petmr_data_path('pet')

def raw_data_path():
    return petmr_data_path('pet')

def _setParameter(hs, set, par, hv):
    h = pystir.cSTIR_setParameter(hs, set, par, hv)
    check_status(h)
    pyiutil.deleteDataHandle(h)
def _set_char_par(handle, set, par, value):
    h = pyiutil.charDataHandle(value)
    _setParameter(handle, set, par, h)
    pyiutil.deleteDataHandle(h)
def _set_int_par(handle, set, par, value):
    h = pyiutil.intDataHandle(value)
    _setParameter(handle, set, par, h)
    pyiutil.deleteDataHandle(h)
def _set_float_par(handle, set, par, value):
    h = pyiutil.floatDataHandle(value)
    _setParameter(handle, set, par, h)
    pyiutil.deleteDataHandle(h)
def _char_par(handle, set, par):
    h = pystir.cSTIR_parameter(handle, set, par)
    check_status(h)
    value = pyiutil.charDataFromHandle(h)
    pyiutil.deleteDataHandle(h)
    return value
def _int_par(handle, set, par):
    h = pystir.cSTIR_parameter(handle, set, par)
    check_status(h)
    value = pyiutil.intDataFromHandle(h)
    pyiutil.deleteDataHandle(h)
    return value
def _float_par(handle, set, par):
    h = pystir.cSTIR_parameter(handle, set, par)
    check_status(h)
    value = pyiutil.floatDataFromHandle(h)
    pyiutil.deleteDataHandle(h)
    return value
def _getParameterHandle(hs, set, par):
    handle = pystir.cSTIR_parameter(hs, set, par)
    check_status(handle)
    return handle

def _tmp_filename():
    return repr(int(1000*time.time()))

class Printer:
    'Redirects STIR printing to files/stdout'
    def __init__(self, info = None, warn = 'stdout', errr = 'stdout'):
        self.info_case = -1
        self.warn_case = -1
        self.errr_case = -1
        if info is None:
            pass
        elif type(info) is not type(' '):
            raise error('wrong info argument for Printer constructor')
        elif info in {'stdout', 'stderr', 'cout', 'cerr'}:
            self.info = pystir.newTextPrinter(info)
            self.info_case = 0
        else:
            self.info = pystir.newTextWriter(info)
            self.info_case = 1
        if warn is None:
            pass
        elif type(warn) is not type(' '):
            raise error('wrong warn argument for Printer constructor')
        elif warn in {'stdout', 'stderr', 'cout', 'cerr'}:
            self.warn = pystir.newTextPrinter(warn)
            self.warn_case = 0
        else:
            self.warn = pystir.newTextWriter(warn)
            self.warn_case = 1
        if errr is None:
            pass
        elif type(errr) is not type(' '):
            raise error('wrong errr argument for Printer constructor')
        elif errr in {'stdout', 'stderr', 'cout', 'cerr'}:
            self.errr = pystir.newTextPrinter(errr)
            self.errr_case = 0
        else:
            self.errr = pystir.newTextWriter(errr)
            self.errr_case = 1
        if self.info_case is not -1:
            pystir.openChannel(0, self.info)
        if self.warn_case is not -1:
            pystir.openChannel(1, self.warn)
        if self.errr_case is not -1:
            pystir.openChannel(2, self.errr)
    def __del__(self):
        if self.info_case is not -1:
            if self.info_case == 0:
                pystir.deleteTextPrinter(self.info)
            else:
                pystir.deleteTextWriter(self.info)
            pystir.closeChannel(0, self.info)
        if self.warn_case is not -1:
            if self.warn_case == 0:
                pystir.deleteTextPrinter(self.warn)
            else:
                pystir.deleteTextWriter(self.warn)
            pystir.closeChannel(1, self.warn)
        if self.errr_case is not -1:
            if self.errr_case == 0:
                pystir.deleteTextPrinter(self.errr)
            else:
                pystir.deleteTextWriter(self.errr)
            pystir.closeChannel(2, self.errr)

class printerTo:
    def __init__(self, dest, channel = -1):
        self.case = -1
        if dest is None:
            return
        elif type(dest) is not type(' '):
            raise error('wrong info argument for Printer constructor')
        elif dest in {'stdout', 'stderr', 'cout', 'cerr'}:
            self.printer = pystir.newTextPrinter(dest)
            self.case = 0
        else:
            self.printer = pystir.newTextWriter(dest)
            self.case = 1
        pystir.openChannel(channel, self.printer)
        self.channel = channel
    def __del__(self):
        if self.case is not -1:
            pystir.closeChannel(self.channel, self.printer)
            if self.case == 0:
                pystir.deleteTextPrinter(self.printer)
            else:
                pystir.deleteTextWriter(self.printer)

class Shape:
    def __init__(self):
        self.handle = None
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def set_origin(self, origin):
        _set_float_par(self.handle, 'Shape', 'x', origin[0])
        _set_float_par(self.handle, 'Shape', 'y', origin[1])
        _set_float_par(self.handle, 'Shape', 'z', origin[2])
    def get_origin(self):
        x = _float_par(self.handle, 'Shape', 'x')
        y = _float_par(self.handle, 'Shape', 'y')
        z = _float_par(self.handle, 'Shape', 'z')
        return (x, y, z)

class EllipsoidalCylinder(Shape):
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

class Voxels:
    def __init__(self, dim, vsize, origin = (0,0,0)):
        self.handle = None
        self.name = 'Voxels'
        self.handle = pystir.cSTIR_voxels3DF\
                      (dim[0], dim[1], dim[2], \
                       vsize[0], vsize[1], vsize[2], \
                       origin[0], origin[1], origin[2])
        check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

class ImageData:
    def __init__(self, arg = None):
        self.handle = None
        if isinstance(arg, str):
            self.handle = pystir.cSTIR_objectFromFile('Image', arg)
            check_status(self.handle)
        elif isinstance(arg, AcquisitionData):
            self.handle = pystir.cSTIR_imageFromAcquisitionData(arg.handle)
            check_status(self.handle)
        elif arg is None:
            self.handle = pystir.cSTIR_newObject('Image')
        else:
            raise error\
                  ('wrong argument ' + repr(arg) + ' for ImageData constructor')
        self.name = 'ImageData'
        self.rimsize = -1
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def initialise\
        (self, arg1, arg2 = 0, arg3 = 0, arg4 = 1, arg5 = 1, arg6 = 1, \
         arg7 = 0, arg8 = 0, arg9 = 0):
        if type(arg1) == type((0,0,0)):
            dim = arg1
            if arg2 == 0:
                vsize = (1, 1, 1)
                origin = (0, 0, 0)
            else:
                vsize = arg2
                if arg3 == 0:
                    origin = (0, 0, 0)
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
        if isinstance(value, numpy.ndarray):
            pystir.cSTIR_setImageData(self.handle, value.ctypes.data)
        else:
            pystir.cSTIR_fillImage(self.handle, value)
        return self
    def clone(self):
        image = ImageData()
        pyiutil.deleteDataHandle(image.handle)
        image.handle = pystir.cSTIR_imageFromImage(self.handle)
        check_status(image.handle)
        return image
    def get_empty_copy(self, value = 1.0):
        image = ImageData()
        pyiutil.deleteDataHandle(image.handle)
        image.handle = pystir.cSTIR_imageFromImage(self.handle)
        check_status(image.handle)
        image.fill(value)
        return image
    def add_shape(self, shape, scale):
        if self.handle is None:
            raise error('cannot add shapes to uninitialised image')
        handle = pystir.cSTIR_addShape(self.handle, shape.handle, scale)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
    def read_from_file(self, filename):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
        self.handle = pystir.cSTIR_objectFromFile('Image', filename)
        check_status(self.handle)
    def diff_from(self, image):
        handle = pystir.cSTIR_imagesDifference\
                 (self.handle, image.handle, self.rimsize)
        check_status(handle)
        diff = pyiutil.doubleDataFromHandle(handle)
        pyiutil.deleteDataHandle(handle)
        return diff
    def as_array(self):
        dim = numpy.ndarray((3,), dtype = numpy.int32)
        pystir.cSTIR_getImageDimensions(self.handle, dim.ctypes.data)
        nz = dim[0]
        ny = dim[1]
        nx = dim[2]
        if nx == 0 or ny == 0 or nz == 0:
            raise error('density data not available')
        array = numpy.ndarray((nz, ny, nx), dtype = numpy.float64)
        pystir.cSTIR_getImageData(self.handle, array.ctypes.data)
        return array
    def show(self):
        if not HAVE_PYLAB:
            print('pylab not found')
            return
        data = self.as_array()
        nz = data.shape[0]
        print('Please enter the number of the slice to view')
        print('(a value outside the range [1 : %d] will stop this loop)' % nz)
        while True:
            s = str(input('slice: '))
            if len(s) < 1:
                break
            z = int(s)
            if z < 1 or z > nz:
                break
            pylab.figure(z)
            pylab.title('image %d' % z)
            pylab.imshow(data[z - 1, :, :])
            print('Close Figure %d window to continue...' % z)
            pylab.show()

class DataProcessor:
    def __init__(self):
        self.handle = None
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def apply(self, image):
        handle = pystir.cSTIR_applyDataProcessor\
                 (self.handle, image.handle)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
    def __del__(self):
        pyiutil.deleteDataHandle(self.handle)

class CylindricFilter(DataProcessor):
    def __init__(self, data_processor = None):
        self.handle = None
        self.name = 'TruncateToCylindricalFOVImageProcessor'
        if data_processor is None:
            self.handle = pystir.cSTIR_newObject(self.name)
        else:
            self.handle = pystir.copyOfObjectHandle(data_processor.handle)
##            self.handle = pystir.cSTIR_copyOfObject(data_processor.handle)
        check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def set_strictly_less_than_radius(self, flag):
        _set_char_par\
            (self.handle, 'TruncateToCylindricalFOVImageProcessor',\
             'strictly_less_than_radius', repr(flag))
    def get_strictly_less_than_radius(self):
        return _int_par\
               (self.handle, 'TruncateToCylindricalFOVImageProcessor',\
                'strictly_less_than_radius') != 0

class RayTracingMatrix:
    def __init__(self):
        self.handle = None
        self.name = 'RayTracingMatrix'
        self.handle = pystir.cSTIR_newObject(self.name)
        check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def set_num_tangential_LORs(self, value):
        _set_int_par(self.handle, self.name, 'num_tangential_LORs', value)
        return self
    def get_num_tangential_LORs(self):
        return _int_par(self.handle, self.name, 'num_tangential_LORs')

class AcquisitionData:
    'Class for PET acquisition data.'
    def __init__(self, src = None):
        ''' Creates new AcquisitionData object from a file or another
            AcquisitionData object;
            src:  file name (Python str) or AcquisitionData object.
        '''
        self.handle = None
        self.name = 'AcquisitionData'
        if src is None:
            return
        if isinstance(src, str):
            self.handle = pystir.cSTIR_objectFromFile('AcquisitionData', src)
        elif isinstance(src, AcquisitionData):
            self.handle = pystir.cSTIR_acquisitionsDataFromTemplate\
                (src.handle)
        else:
            raise error('wrong source in AcquisitionData constructor')
        check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def create_empty_image(self, value = 0):
        ''' Creates ImageData object containing PET image of dimensions
            and voxel sizes compatible with the scanner geometry stored
            in this AcquisitionData object and assigns a given value
            to all voxels;
            value:  a Python float.
        '''
        image = ImageData()
        pyiutil.deleteDataHandle(image.handle)
        image.handle = pystir.cSTIR_imageFromAcquisitionData(self.handle)
        check_status(image.handle)
        image.fill(value)
        return image
    def as_array(self):
        ''' Returns a copy of acquisition data stored in this object as a
            NumPy ndarray of 3 dimensions (in default C ordering of data):
            - number of sinograms
            - number of views
            - number of tangential positions.
        '''
        dim = numpy.ndarray((3,), dtype = numpy.int32)
        handle = pystir.cSTIR_getAcquisitionsDimensions\
            (self.handle, dim.ctypes.data)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
        nx = dim[0]
        ny = dim[1]
        nz = dim[2]
        if nx == 0 or ny == 0 or nz == 0:
            raise error('density data not available')
        array = numpy.ndarray((nz, ny, nx), dtype = numpy.float64)
        handle = pystir.cSTIR_getAcquisitionsData\
            (self.handle, array.ctypes.data)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
        return array
    def fill(self, value):
        ''' Fills the object with values;
            value:  either NumPy ndarray or another AcquisitionData object
                    or Python float.
        '''
        if self.handle is None:
            raise error('AcquisitionData object not initialized')
        if isinstance(value, numpy.ndarray):
            pystir.cSTIR_setAcquisitionsData(self.handle, value.ctypes.data)
        elif isinstance(value, AcquisitionData):
            pystir.cSTIR_fillAcquisitionsDataFromAcquisitionsData\
                (self.handle, value.handle)
        elif isinstance(value, float):
            pystir.cSTIR_fillAcquisitionsData(self.handle, value)
        else:
            raise error('wrong fill value')
        return self
    def clone(self):
        ''' Returns a true copy of this object (not Python handle).
        '''
        ad = AcquisitionData(self)
        ad.fill(self)
        return ad
    def get_empty_copy(self, value = 0):
        ''' Returns a copy of this object filled with a given value;
            value:  a Python float.
        '''
        ad = AcquisitionData(self)
        ad.fill(value)
        return ad

class AcquisitionModel:
    ''' Class for a PET acquisition model that relates an image x to the
        acquisition data y as
        (F)    y = [1/n](G x + [a]) + [b]
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
        where G' is the transpose of G, is referred to as backward projection.
    '''
    def __init__(self):
        self.handle = None
        self.name = 'AcquisitionModel'
    def set_up(self, acq_templ, img_templ):
        ''' Prepares this object for performing forward and backward
            projections;
            acq_templ:  an AcquisitionData object used as a template for
                        creating an AcquisitionData object to store forward
                        projection;
            img_templ:  an ImageData object used as a template for creating an
                        ImageData object to store backward projection.
        '''
        handle = pystir.cSTIR_setupAcquisitionModel\
            (self.handle, acq_templ.handle, img_templ.handle)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
    def set_additive_term(self, at):
        ''' Sets the additive term a in (F);
            at:  an AcquisitionData object containing a.
        '''
        _setParameter\
            (self.handle, 'AcquisitionModel', 'additive_term', at.handle)
    def set_background_term(self, bt):
        ''' Sets the additive term b in (F);
            bt:  an AcquisitionData object containing b.
        '''
        _setParameter\
            (self.handle, 'AcquisitionModel', 'background_term', bt.handle)
    def set_normalisation(self, bin_eff):
        ''' Sets the normalization n in (F);
            bin_eff:  an AcquisitionData object containing bin efficiencies
                      (the inverse of n).
        '''
        _setParameter\
            (self.handle, 'AcquisitionModel', 'normalisation', bin_eff.handle)
    def forward(self, image, filename = ''):
        ''' Returns the forward projection of x given by (F);
            image   :  an ImageData object containing x;
            filename:  an optional name of the file to store projection data;
                       if not present, projection data is stored in memory
                       (not recommended as it can be huge).
        '''
        ad = AcquisitionData()
        ad.handle = pystir.cSTIR_acquisitionModelFwd\
            (self.handle, image.handle, filename)
        check_status(ad.handle)
        return ad;
    def backward(self, ad):
        ''' Returns the backward projection of y giben by (B);
            ad:  an AcquisitionData object containing y.
        '''
        image = ImageData()
        image.handle = pystir.cSTIR_acquisitionModelBwd\
            (self.handle, ad.handle)
        check_status(image.handle)
        return image

class AcquisitionModelUsingMatrix(AcquisitionModel):
    ''' Class for a PET acquisition model that uses (implicitly) a sparse
        matrix for G in (F).
    '''
    def __init__(self, matrix = None):
        ''' Creates an AcquisitionModelUsingMatrix object, optionally setting
            the ray tracing matrix to be used for projecting;
            matrix:  a RayTracingMatrix object to represent G in (F).
        '''
        self.handle = None
        self.name = 'AcqModUsingMatrix'
        self.handle = pystir.cSTIR_newObject(self.name)
        check_status(self.handle)
        if matrix is not None:
            _setParameter(self.handle, self.name, 'matrix', matrix.handle)
    def set_matrix(self, matrix):
        ''' Sets the ray tracing matrix to be used for projecting;
            matrix:  a RayTracingMatrix object to represent G in (F).
        '''
        _setParameter(self.handle, self.name, 'matrix', matrix.handle)
    def get_matrix(self):
        ''' Returns the ray tracing matrix used for projecting;
            matrix:  a RayTracingMatrix object representing G in (F).
        '''
        matrix = RayTracingMatrix()
        matrix.handle = pystir.cSTIR_parameter(self.handle, self.name, 'matrix')
        check_status(matrix.handle)
        return matrix

class Prior:
    def __init__(self):
        self.handle = None
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def set_penalisation_factor(self, value):
        _set_float_par\
            (self.handle, 'GeneralisedPrior', 'penalisation_factor', value)
        return self
    def get_penalisation_factor(self):
        return _float_par\
            (self.handle, 'GeneralisedPrior', 'penalisation_factor')
    def get_gradient(self, image):
        grad = ImageData()
        pyiutil.deleteDataHandle(grad.handle)
        grad.handle = pystir.cSTIR_priorGradient(self.handle, image.handle)
        check_status(grad.handle)
        return grad
##    def set_up(self):
##        handle = pystir.cSTIR_setupObject('GeneralisedPrior', self.handle)
##        check_status(handle)
##        pyiutil.deleteDataHandle(handle)

class QuadraticPrior(Prior):
    def __init__(self):
        self.handle = None
        self.name = 'QuadraticPrior'
        self.handle = pystir.cSTIR_newObject(self.name)
        check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)

class ObjectiveFunction:
    def __init__(self):
        self.handle = None
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def set_prior(self, prior):
        _setParameter(self.handle, 'GeneralisedObjectiveFunction',\
            'prior', prior.handle)
        self.prior = prior
    def get_prior(self):
        prior = Prior()
        prior.handle = pystir.cSTIR_parameter\
            (self.handle, 'GeneralisedObjectiveFunction', 'prior')
        check_status(prior.handle)
        return prior
    def set_num_subsets(self, n):
        _set_int_par\
            (self.handle, 'GeneralisedObjectiveFunction', 'num_subsets', n)
    def set_up(self, image):
        handle = pystir.cSTIR_setupObjectiveFunction(self.handle, image.handle)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
    def value(self, image):
        handle = pystir.cSTIR_objectiveFunctionValue(self.handle, image.handle)
        check_status(handle)
        return pyiutil.floatDataFromHandle(handle)
    def gradient(self, image, subset):
        grad = ImageData()
        pyiutil.deleteDataHandle(grad.handle)
        grad.handle = pystir.cSTIR_objectiveFunctionGradient\
            (self.handle, image.handle, subset)
        check_status(grad.handle)
        return grad

class PoissonLogLh_LinModMean(ObjectiveFunction):
    def __init__(self):
        self.handle = None
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def set_sensitivity_filename(self, name):
        _set_char_par\
            (self.handle, 'PoissonLogLikelihoodWithLinearModelForMean',\
             'sensitivity_filename', name)
    def set_use_subset_sensitivities(self, flag):
        _set_char_par\
            (self.handle, 'PoissonLogLikelihoodWithLinearModelForMean',\
             'use_subset_sensitivities', repr(flag))
    def set_recompute_sensitivity(self, flag):
        _set_char_par\
            (self.handle, 'PoissonLogLikelihoodWithLinearModelForMean',\
             'recompute_sensitivity', repr(flag))
    def get_subset_sensitivity(self, subset):
        ss = ImageData()
        pyiutil.deleteDataHandle(ss.handle)
        ss.handle = pystir.cSTIR_subsetSensitivity(self.handle, subset)
        check_status(ss.handle)
        return ss
##    def get_gradient_not_divided(self, image, subset):
    def get_gradient_plus_sensitivity_no_penalty(self, image, subset):
        grad = ImageData()
        pyiutil.deleteDataHandle(grad.handle)
        grad.handle = pystir.cSTIR_objectiveFunctionGradientNotDivided\
            (self.handle, image.handle, subset)
        check_status(grad.handle)
        return grad

class PoissonLogLh_LinModMean_AcqMod(PoissonLogLh_LinModMean):
    def __init__(self, obj_fun = None):
        self.handle = None
        self.name = 'PoissonLogLikelihoodWithLinearModelForMeanAndProjData'
        if obj_fun is None:
            self.handle = pystir.cSTIR_newObject(self.name)
        else:
            self.handle = pyiutil.copyOfObjectHandle(obj_fun.handle)
##            self.handle = pystir.cSTIR_copyOfObject(obj_fun.handle)
        check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def set_input_filename(self, name):
        _set_char_par\
            (self.handle, self.name, 'input_filename', name)
    def set_zero_seg0_end_planes(self, flag):
        _set_char_par\
            (self.handle, self.name, 'zero_seg0_end_planes', repr(flag))
    def set_max_segment_num_to_process(self, n):
        _set_int_par(self.handle, self.name, 'max_segment_num_to_process', n)
    def set_acquisition_model(self, am):
        _setParameter\
            (self.handle, self.name, 'acquisition_model', am.handle)
    def get_acquisition_model(self):
        am = AcquisitionModelUsingMatrix()
        if am.handle is not None:
            pyiutil.deleteDataHandle(am.handle)
        am.handle = pystir.cSTIR_parameter\
            (self.handle, self.name, 'acquisition_model')
        check_status(am.handle)
        return am
    def set_acquisition_data(self, ad):
        _setParameter\
            (self.handle, self.name, 'proj_data_sptr', ad.handle)

class Reconstruction:
    def __init__(self):
        self.handle = None
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def set_output_filename_prefix(self, prefix):
        _set_char_par\
            (self.handle, 'Reconstruction', 'output_filename_prefix', prefix)

class IterativeReconstruction(Reconstruction):
    def __init__(self):
        self.handle = None
        self.input = None
        self.image = None
        self.subset = 0
    def __del__(self):
        if self.handle is not None:
            pystir.cSTIR_deleteObject(self.handle)
    def set_num_subsets(self, n):
        _set_int_par\
            (self.handle, 'IterativeReconstruction', 'num_subsets', n)
    def get_num_subsets(self):
        return _int_par\
            (self.handle, 'IterativeReconstruction', 'num_subsets')
    def set_start_subset_num(self, n):
        _set_int_par\
            (self.handle, 'IterativeReconstruction', 'start_subset_num', n)
    def get_start_subset_num(self):
        return _int_par\
            (self.handle, 'IterativeReconstruction', 'start_subset_num')
    def set_num_subiterations(self, n):
        _set_int_par\
            (self.handle, 'IterativeReconstruction', 'num_subiterations', n)
    def get_num_subiterations(self):
        return _int_par\
            (self.handle, 'IterativeReconstruction', 'num_subiterations')
    def set_start_subiteration_num(self, n):
        _set_int_par\
            (self.handle, 'IterativeReconstruction', 'start_subiteration_num', n)
    def get_start_subiteration_num(self):
        return _int_par\
            (self.handle, 'IterativeReconstruction', 'start_subiteration_num')
    def set_subiteration_num(self, iter):
        _set_int_par\
            (self.handle, 'IterativeReconstruction', 'subiteration_num', iter)
    def get_subiteration_num(self):
        return _int_par\
            (self.handle, 'IterativeReconstruction', 'subiteration_num')
    def set_save_interval(self, n):
        _set_int_par\
            (self.handle, 'IterativeReconstruction', 'save_interval', n)
    def set_inter_iteration_filter_interval(self, n):
        _set_int_par\
            (self.handle, 'IterativeReconstruction',\
             'inter_iteration_filter_interval', n)
    def set_objective_function(self, obj):
        _setParameter\
            (self.handle, 'IterativeReconstruction',\
             'objective_function', obj.handle)
    def get_objective_function(self):
        obj_fun = ObjectiveFunction()
        obj_fun.handle = pystir.cSTIR_parameter\
            (self.handle, 'IterativeReconstruction', 'objective_function')
        check_status(obj_fun.handle)
        return obj_fun
    def set_inter_iteration_filter(self, f):
        pystir.cSTIR_setParameter\
            (self.handle, 'IterativeReconstruction',\
             'inter_iteration_filter_type', f.handle)
    def get_inter_iteration_filter(self):
        filter = DataProcessor()
        filter.handle = pystir.cSTIR_parameter\
            (self.handle, 'IterativeReconstruction',\
             'inter_iteration_filter_type')
        check_status(filter.handle)
        return filter
    def set_up(self, image):
        handle = pystir.cSTIR_setupReconstruction(self.handle, image.handle)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
    def set_input(self, input_data):
        self.input = input_data
    def set_current_estimate(self, image):
        self.image = image
    def get_current_estimate(self):
        return self.image
    def update_current_estimate(self):
        if self.image is None:
            raise error('current estimate not set')
        handle = pystir.cSTIR_updateReconstruction(self.handle, self.image.handle)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
    def set_current_subset(self, subset):
        self.subset = subset
    def get_subset_sensitivity(self):
        obj_fun = self.get_objective_function()
        obj_fun.get_subset_sensitivity(self.subset)
    def reconstruct(self, image):
        handle = pystir.cSTIR_runReconstruction(self.handle, image.handle)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
    def process(self):
        if self.image is None:
            raise error('current estimate not set')
        handle = pystir.cSTIR_runReconstruction(self.handle, self.image.handle)
        check_status(handle)
        pyiutil.deleteDataHandle(handle)
    def update(self, image):
        self.set_current_estimate(image)
        self.update_current_estimate()
        return self.get_current_estimate()

class OSMAPOSLReconstruction(IterativeReconstruction):
    def __init__(self, filename = ''):
        self.handle = None
        self.name = 'OSMAPOSL'
        self.handle = pystir.cSTIR_objectFromFile\
            ('OSMAPOSLReconstruction', filename)
        check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pyiutil.deleteDataHandle(self.handle)
    def set_MAP_model(self, model):
        _set_char_par\
            (self.handle, self.name, 'MAP_model', model)
    def get_objective_function(self):
        obj_fun = PoissonLogLh_LinModMean()
        obj_fun.handle = pystir.cSTIR_parameter\
            (self.handle, self.name, 'objective_function')
        check_status(obj_fun.handle)
        return obj_fun

class OSSPSReconstruction(IterativeReconstruction):
    def __init__(self, filename = ''):
        self.handle = None
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
