import numpy
import os
import pystir
import time

INFO_CHANNEL = 0
WARNING_CHANNEL = 1
ERROR_CHANNEL = 2
ALL_CHANNELS = -1

class error(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def _check_status(handle):
    if pystir.executionStatus(handle) != 0:
        msg = pystir.executionError(handle)
        file = pystir.executionErrorFile(handle)
        line = pystir.executionErrorLine(handle)
        errorMsg = \
            repr(msg) + ' exception thrown at line ' + \
            repr(line) + ' of ' + file
        raise error(errorMsg)

def _setParameter(hs, set, par, hv):
    h = pystir.cSTIR_setParameter(hs, set, par, hv)
    _check_status(h)
    pystir.deleteDataHandle(h)
def _set_char_par(handle, set, par, value):
    h = pystir.charDataHandle(value)
    _setParameter(handle, set, par, h)
    pystir.deleteDataHandle(h)
def _set_int_par(handle, set, par, value):
    h = pystir.intDataHandle(value)
    _setParameter(handle, set, par, h)
    pystir.deleteDataHandle(h)
def _set_float_par(handle, set, par, value):
    h = pystir.floatDataHandle(value)
    _setParameter(handle, set, par, h)
    pystir.deleteDataHandle(h)
def _char_par(handle, set, par):
    h = pystir.cSTIR_parameter(handle, set, par)
    _check_status(h)
    value = pystir.charDataFromHandle(h)
    pystir.deleteDataHandle(h)
    return value
def _int_par(handle, set, par):
    h = pystir.cSTIR_parameter(handle, set, par)
    _check_status(h)
    value = pystir.intDataFromHandle(h)
    pystir.deleteDataHandle(h)
    return value
def _float_par(handle, set, par):
    h = pystir.cSTIR_parameter(handle, set, par)
    _check_status(h)
    value = pystir.floatDataFromHandle(h)
    pystir.deleteDataHandle(h)
    return value
def _getParameterHandle(hs, set, par):
    handle = pystir.cSTIR_parameter(hs, set, par)
    _check_status(handle)
    return handle

def _tmp_filename():
    return repr(int(1000*time.time()))

class Printer:
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
            pystir.deleteDataHandle(self.handle)
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
        _check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pystir.deleteDataHandle(self.handle)
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
        _check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pystir.deleteDataHandle(self.handle)

class PETImage:
    def __init__(self, arg = None):
        self.handle = None
        self.voxels = None
        if isinstance(arg, str):
            self.handle = pystir.cSTIR_objectFromFile('Image', arg)
            _check_status(self.handle)
        elif arg is None:
            self.handle = pystir.cSTIR_newObject('Image')
        else:
            raise error\
                  ('wrong argument ' + repr(arg) + ' for Image constructor')
        self.name = 'Image'
        self.rimsize = -1
    def __del__(self):
        if self.handle is not None:
            pystir.deleteDataHandle(self.handle)
        if self.voxels is not None:
            pystir.deleteDataHandle(self.voxels)
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
            pystir.deleteDataHandle(self.handle)
        if self.voxels is not None:
            pystir.deleteDataHandle(self.voxels)
        self.handle = None
        self.voxels = pystir.cSTIR_voxels3DF\
                      (dim[0], dim[1], dim[2], \
                       vsize[0], vsize[1], vsize[2], \
                       origin[0], origin[1], origin[2])
        _check_status(self.voxels)
        self.handle = pystir.cSTIR_imageFromVoxels(self.voxels)
        _check_status(self.handle)
    def fill(self, value):
        if isinstance(value, numpy.ndarray):
            pystir.cSTIR_setImageData(self.handle, value.ctypes.data)
        else:
            pystir.cSTIR_fillImage(self.handle, value)
        return self
    def clone(self):
        image = PETImage()
        image.handle = pystir.cSTIR_imageFromImage(self.handle)
        _check_status(image.handle)
        return image
    def get_empty_copy(self, value = 1.0):
        image = PETImage()
        image.handle = pystir.cSTIR_imageFromImage(self.handle)
        _check_status(image.handle)
        image.fill(value)
        return image
    def add_shape(self, shape, scale):
        if self.handle is None:
            raise error('cannot add shapes to uninitialised image')
        handle = pystir.cSTIR_addShape\
                 (self.handle, self.voxels, shape.handle, scale)
        _check_status(handle)
        pystir.deleteDataHandle(handle)
    def read_from_file(self, filename):
        if self.handle is not None:
            pystir.deleteDataHandle(self.handle)
        if self.voxels is not None:
            pystir.deleteDataHandle(self.voxels)
        self.handle = pystir.cSTIR_objectFromFile('Image', filename)
        _check_status(self.handle)
    def diff_from(self, image):
        handle = pystir.cSTIR_imagesDifference\
                 (self.handle, image.handle, self.rimsize)
        _check_status(handle)
        diff = pystir.doubleDataFromHandle(handle)
        pystir.deleteDataHandle(handle)
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

class DataProcessor:
    def __init__(self):
        self.handle = None
    def __del__(self):
        if self.handle is not None:
            pystir.deleteDataHandle(self.handle)
    def apply(self, image):
        handle = pystir.cSTIR_applyDataProcessor\
                 (self.handle, image.handle)
        _check_status(handle)
        pystir.deleteDataHandle(handle)
    def __del__(self):
        pystir.deleteDataHandle(self.handle)

class CylindricFilter(DataProcessor):
    def __init__(self, data_processor = None):
        self.handle = None
        self.name = 'TruncateToCylindricalFOVImageProcessor'
        if data_processor is None:
            self.handle = pystir.cSTIR_newObject(self.name)
        else:
            self.handle = pystir.copyOfObjectHandle(data_processor.handle)
##            self.handle = pystir.cSTIR_copyOfObject(data_processor.handle)
        _check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pystir.deleteDataHandle(self.handle)
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
        _check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pystir.deleteDataHandle(self.handle)
    def set_num_tangential_LORs(self, value):
        _set_int_par(self.handle, self.name, 'num_tangential_LORs', value)
    def get_num_tangential_LORs(self):
        return _int_par(self.handle, self.name, 'num_tangential_LORs')

class PETAcquisitionData:
    def __init__(self, src = None):
        self.handle = None
        self.name = 'AcquisitionData'
        if src is None:
            return
        if isinstance(src, str):
            self.handle = pystir.cSTIR_objectFromFile('AcquisitionData', src)
        elif isinstance(src, PETAcquisitionData):
            self.handle = pystir.cSTIR_acquisitionsDataFromTemplate\
                (src.handle)
        else:
            raise error('wrong source in AcquisitionData constructor')
        _check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pystir.deleteDataHandle(self.handle)
    def as_array(self):
        dim = numpy.ndarray((3,), dtype = numpy.int32)
        handle = pystir.cSTIR_getAcquisitionsDimensions\
            (self.handle, dim.ctypes.data)
        _check_status(handle)
        pystir.deleteDataHandle(handle)
        nx = dim[0]
        ny = dim[1]
        nz = dim[2]
        if nx == 0 or ny == 0 or nz == 0:
            raise error('density data not available')
        array = numpy.ndarray((nz, ny, nx), dtype = numpy.float64)
        handle = pystir.cSTIR_getAcquisitionsData\
            (self.handle, array.ctypes.data)
        _check_status(handle)
        pystir.deleteDataHandle(handle)
        return array
    def fill(self, value):
        if isinstance(value, numpy.ndarray):
            pystir.cSTIR_setAcquisitionsData(self.handle, value.ctypes.data)
        elif isinstance(value, PETAcquisitionData):
            pystir.cSTIR_fillAcquisitionsDataFromAcquisitionsData\
                (self.handle, value)
        else : # should check on double really
            pystir.cSTIR_fillAcquisitionsData(self.handle, value)
        return self
    def clone(self):
        ad = PETAcquisitionData()
        ad.handle = pystir.cSTIR_acquisitionsDataFromTemplate(self.handle)
        _check_status(ad.handle)
        ad.fill(self)
        return ad
    def get_empty_copy(self, value = 0):
        ad = PETAcquisitionData()
        ad.handle = pystir.cSTIR_acquisitionsDataFromTemplate(self.handle)
        _check_status(ad.handle)
        ad.fill(value)
        return ad

class PETAcquisitionModel:
    def __init__(self):
        self.handle = None
        self.name = 'AcquisitionModel'
    def set_up(self, templ, image):
        handle = pystir.cSTIR_setupAcquisitionModel\
            (self.handle, templ.handle, image.handle)
        _check_status(handle)
    def set_additive_term(self, at):
        _setParameter\
            (self.handle, 'AcquisitionModel', 'additive_term', at.handle)
    def set_background_term(self, at):
        _setParameter\
            (self.handle, 'AcquisitionModel', 'background_term', at.handle)
    def set_normalisation(self, norm):
        _setParameter\
            (self.handle, 'AcquisitionModel', 'normalisation', norm.handle)
    def forward(self, image, filename = ''):
        ad = PETAcquisitionData()
        ad.handle = pystir.cSTIR_acquisitionModelFwd\
            (self.handle, image.handle, filename)
        _check_status(ad.handle)
        return ad;
    def backward(self, ad):
        image = PETImage()
        image.handle = pystir.cSTIR_acquisitionModelBwd\
            (self.handle, ad.handle)
        _check_status(image.handle)
        return image

class PETAcquisitionModelUsingMatrix(PETAcquisitionModel):
    def __init__(self):
        self.handle = None
        self.name = 'AcqModUsingMatrix'
        self.handle = pystir.cSTIR_newObject(self.name)
        _check_status(self.handle)
    def set_matrix(self, matrix):
        _setParameter(self.handle, self.name, 'matrix', matrix.handle)
    def get_matrix(self):
        matrix = RayTracingMatrix()
        matrix.handle = pystir.cSTIR_parameter\
            (self.handle, self.name, 'matrix')
        _check_status(matrix.handle)
        return matrix

class Prior:
    def __init__(self):
        self.handle = None
    def __del__(self):
        if self.handle is not None:
            pystir.deleteDataHandle(self.handle)
    def set_penalisation_factor(self, value):
        _set_float_par\
            (self.handle, 'GeneralisedPrior', 'penalisation_factor', value)
    def get_penalisation_factor(self):
        return _float_par\
        (self.handle, 'GeneralisedPrior', 'penalisation_factor')
    def set_up(self):
        handle = pystir.cSTIR_setupObject('GeneralisedPrior', self.handle)
        _check_status(handle)
        pystir.deleteDataHandle(handle)

class QuadraticPrior(Prior):
    def __init__(self):
        self.handle = None
        self.name = 'QuadraticPrior'
        self.handle = pystir.cSTIR_newObject(self.name)
        _check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pystir.deleteDataHandle(self.handle)

class ObjectiveFunction:
    def __init__(self):
        self.handle = None
    def __del__(self):
        if self.handle is not None:
            pystir.deleteDataHandle(self.handle)
    def set_prior(self, prior):
        _setParameter(self.handle, 'GeneralisedObjectiveFunction',\
            'prior', prior.handle)
        self.prior = prior
    def get_prior(self):
        prior = Prior()
        prior.handle = pystir.cSTIR_parameter\
            (self.handle, 'GeneralisedObjectiveFunction', 'prior')
        _check_status(prior.handle)
        return prior
    def set_up(self):
        handle = pystir.cSTIR_setupObject\
            ('GeneralisedObjectiveFunction', self.handle)
        _check_status(handle)
        pystir.deleteDataHandle(handle)
    def value(self, image):
        handle = pystir.cSTIR_objectiveFunctionValue(self.handle, image.handle)
        _check_status(handle)
        return pystir.floatDataFromHandle(handle)
    def gradient(self, image, subset):
        grad = PETImage()
        grad.handle = pystir.cSTIR_objectiveFunctionGradient\
            (self.handle, image.handle, subset)
        _check_status(grad.handle)
        return grad

class PoissonLogLh_LinModMean(ObjectiveFunction):
    def __init__(self):
        self.handle = None
    def __del__(self):
        if self.handle is not None:
            pystir.deleteDataHandle(self.handle)
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

class PoissonLogLh_LinModMean_AcqModData(PoissonLogLh_LinModMean):
    def __init__(self, obj_fun = None):
        self.handle = None
        self.name = 'PoissonLogLikelihoodWithLinearModelForMeanAndProjData'
        if obj_fun is None:
            self.handle = pystir.cSTIR_newObject(self.name)
        else:
            self.handle = pystir.copyOfObjectHandle(obj_fun.handle)
##            self.handle = pystir.cSTIR_copyOfObject(obj_fun.handle)
        _check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pystir.deleteDataHandle(self.handle)
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
        am = PETAcquisitionModelUsingMatrix()
        if am.handle is not None:
            pystir.deleteDataHandle(am.handle)
        am.handle = pystir.cSTIR_parameter\
            (self.handle, self.name, 'acquisition_model')
        _check_status(am.handle)
        return am
    def set_acquisition_data(self, am):
        _setParameter\
            (self.handle, self.name, 'proj_data_sptr', am.handle)

class Reconstruction:
    def __init__(self):
        self.handle = None
    def __del__(self):
        if self.handle is not None:
            pystir.deleteDataHandle(self.handle)
    def set_output_filename_prefix(self, prefix):
        _set_char_par\
            (self.handle, 'Reconstruction', 'output_filename_prefix', prefix)

class IterativeReconstruction(Reconstruction):
    def __init__(self):
        self.handle = None
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
        _check_status(obj_fun.handle)
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
        _check_status(filter.handle)
        return filter
    def set_up(self, image):
        handle = pystir.cSTIR_setupReconstruction(self.handle, image.handle)
        _check_status(handle)
        pystir.deleteDataHandle(handle)
    def reconstruct(self, image):
        handle = pystir.cSTIR_runReconstruction(self.handle, image.handle)
        _check_status(handle)
        pystir.deleteDataHandle(handle)
    def update(self, image):
        handle = pystir.cSTIR_updateReconstruction(self.handle, image.handle)
        _check_status(handle)
        pystir.deleteDataHandle(handle)

class OSMAPOSLReconstruction(IterativeReconstruction):
    def __init__(self, filename = ''):
        self.handle = None
        self.name = 'OSMAPOSL'
        self.handle = pystir.cSTIR_objectFromFile\
            ('OSMAPOSLReconstruction', filename)
        _check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pystir.deleteDataHandle(self.handle)
    def set_MAP_model(self, model):
        _set_char_par\
            (self.handle, self.name, 'MAP_model', model)
    def get_objective_function(self):
        obj_fun = PoissonLogLh_LinModMean()
        obj_fun.handle = pystir.cSTIR_parameter\
            (self.handle, self.name, 'objective_function')
        _check_status(obj_fun.handle)
        return obj_fun

class OSSPSReconstruction(IterativeReconstruction):
    def __init__(self, filename = ''):
        self.handle = None
        self.name = 'OSSPS'
        self.handle = pystir.cSTIR_objectFromFile\
                      ('OSSPSReconstruction', filename)
        _check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pystir.deleteDataHandle(self.handle)
    def set_relaxation_parameter(self, value):
        _set_float_par\
            (self.handle, self.name, 'relaxation_parameter', value)
