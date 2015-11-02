import numpy
import pystir

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
        errorMsg =\
            msg + ' exception thrown at line ' + repr(line) + ' of ' + file
        raise error(errorMsg)

def _setParameter(hs, set, par, hv):
    h = pystir.cSTIR_setParameter(hs, set, par, hv)
    _check_status(h)
    pystir.deleteDataHandle(h)

def set_int_par(handle, set, par, value):
    h = pystir.intDataHandle(value)
    _setParameter(handle, set, par, h)
    #pystir.cSTIR_setParameter(handle, set, par, h)
    pystir.deleteDataHandle(h)
##    if status != 0:
##        errorMsg = 'parameter ' + par + ' not found'
##        raise error(errorMsg)
def get_int_par(handle, set, par):
    h = pystir.cSTIR_parameter(handle, set, par)
    _check_status(h)
    value = pystir.intDataFromHandle(h)
    pystir.deleteDataHandle(h)
    return value

class printerTo:
    def __init__(self, dest):
        handle = pystir.newTextWriterHandle()
        #if dest.casefold() in {'stdout', 'stderr', 'cout', 'cerr'}:
        if dest in {'stdout', 'stderr', 'cout', 'cerr'}:
            self.printer = pystir.newTextPrinter(dest)
            self.type = 0
        else:
            self.printer = pystir.newTextWriter(dest)
            self.type = 1
        pystir.setWriter(handle, self.printer)
        pystir.deleteTextWriterHandle(handle)
    def __del__(self):
        if self.type == 0:
            pystir.deleteTextPrinter(self.printer)
        else:
            pystir.deleteTextWriter(self.printer)
        pystir.resetWriter()

class Image3DFloat:
    def __init__(self):
        self.handle = pystir.newDataHandle()
        self.name = 'Image3DF'
        self.rimsize = -1
    def __del__(self):
        pystir.deleteDataHandle(self.handle)
    def read_from_file(self, filename):
        pystir.deleteDataHandle(self.handle)
        self.handle = pystir.cSTIR_imageFromFile(filename)
        _check_status(self.handle)
    def diff_from(self, image):
        handle = pystir.cSTIR_imagesDifference\
                 (self.handle, image.handle, self.rimsize)
        _check_status(handle)
        diff = pystir.doubleDataFromHandle(handle)
        pystir.deleteDataHandle(handle)
        return diff

def imageData(image):
    dim = numpy.ndarray((3,), dtype = numpy.int32)
    pystir.cSTIR_getImageDimensions(image.handle, dim.ctypes.data)
    nz = dim[0]
    ny = dim[1]
    nx = dim[2]
    data = numpy.ndarray((nz, ny, nx), dtype = numpy.float64)
    pystir.cSTIR_getImageData(image.handle, data.ctypes.data)
    return data

class TruncateToCylindricalFOVImageProcessor:
    def __init__(self):
        self.name = 'TruncateToCylindricalFOVImageProcessor'
        self.handle = pystir.cSTIR_newObject(self.name);
    def __del__(self):
        pystir.deleteDataHandle(self.handle)

class ProjectorsUsingMatrix:
    def __init__(self):
        self.name = 'ProjectorsUsingMatrix'
        self.handle = pystir.cSTIR_newObject(self.name);
    def __del__(self):
        pystir.deleteDataHandle(self.handle)
    def set_matrix_type(self, value):
        #pystir.cSTIR_setParameter\
        _setParameter\
            (self.handle, 'ProjectorsUsingMatrix', 'matrix_type', value.handle)

class RayTracingMatrix:
    def __init__(self):
        self.name = 'RayTracingMatrix'
        self.handle = pystir.cSTIR_newObject(self.name);
    def __del__(self):
        pystir.deleteDataHandle(self.handle)
    def set_num_tangential_LORs(self, value):
        h = pystir.intDataHandle(value)
        #pystir.cSTIR_setParameter\
        _setParameter\
            (self.handle, 'RayTracingMatrix', 'num_tangential_LORs', h)
        pystir.deleteDataHandle(h)
    def get_num_tangential_LORs(self):
        h = pystir.cSTIR_parameter\
            (self.handle, 'RayTracingMatrix', 'num_tangential_LORs')
        value = pystir.intDataFromHandle(h)
        pystir.deleteDataHandle(h)
        return value

class GeneralisedPrior:
    def set_penalisation_factor(self, value):
        h = pystir.floatDataHandle(value)
        #pystir.cSTIR_setParameter\
        _setParameter\
            (self.handle, 'GeneralisedPrior', 'penalisation_factor', h)
        pystir.deleteDataHandle(h)
    def get_penalisation_factor(self):
        h = pystir.cSTIR_parameter\
            (self.handle, 'GeneralisedPrior', 'penalisation_factor')
        value = pystir.floatDataFromHandle(h)
        pystir.deleteDataHandle(h)
        return value
    def set_up(self):
        handle = pystir.cSTIR_setupObject("prior", self.handle)
        _check_status(handle)
        pystir.deleteDataHandle(handle)

class QuadraticPrior(GeneralisedPrior):
    def __init__(self):
        self.name = 'QuadraticPrior'
        self.handle = pystir.cSTIR_newObject(self.name);
    def __del__(self):
        pystir.deleteDataHandle(self.handle)

class GeneralisedObjectiveFunction:
    def __init__(self):
        self.handle = pystir.newDataHandle()
        self.name = 'GeneralisedObjectiveFunction'
    def __del__(self):
        pystir.deleteDataHandle(self.handle)
    def set_prior(self, prior):
        #pystir.cSTIR_setParameter\
        _setParameter\
            (self.handle, 'GeneralisedObjectiveFunction', 'prior', prior.handle)
    def set_up(self):
        handle = pystir.cSTIR_setupObject("objective_function", self.handle)
        _check_status(handle)
        pystir.deleteDataHandle(handle)

class PoissonLogLikelihoodWithLinearModelForMean3DFloat\
      (GeneralisedObjectiveFunction):
    def __init__(self):
        self.handle = pystir.newDataHandle()
        self.name = 'PoissonLogLikelihoodWithLinearModelForMean'
    def __del__(self):
        pystir.deleteDataHandle(self.handle)
    def set_sensitivity_filename(self, name):
        handle = pystir.charDataHandle(name)
        #pystir.cSTIR_setParameter\
        _setParameter\
            (self.handle, 'PoissonLogLikelihoodWithLinearModelForMean',\
             'sensitivity_filename', handle)
        pystir.deleteDataHandle(handle)
    def set_use_subset_sensitivities(self, flag):
        handle = pystir.charDataHandle(repr(flag))
        #pystir.cSTIR_setParameter\
        _setParameter\
            (self.handle, 'PoissonLogLikelihoodWithLinearModelForMean',\
             'use_subset_sensitivities', handle)
        pystir.deleteDataHandle(handle)
    def set_recompute_sensitivity(self, value):
        handle = pystir.charDataHandle(repr(value))
        #pystir.cSTIR_setParameter\
        _setParameter\
            (self.handle, 'PoissonLogLikelihoodWithLinearModelForMean',\
             'recompute_sensitivity', handle)
        pystir.deleteDataHandle(handle)

class PoissonLogLikelihoodWithLinearModelForMeanAndProjData3DFloat\
      (PoissonLogLikelihoodWithLinearModelForMean3DFloat):
    def __init__(self):
        self.name = 'PoissonLogLikelihoodWithLinearModelForMeanAndProjData'
        self.handle = pystir.cSTIR_newObject(self.name);
    def __del__(self):
        pystir.deleteDataHandle(self.handle)
    def set_input_filename(self, name):
        handle = pystir.charDataHandle(name)
        #pystir.cSTIR_setParameter\
        _setParameter\
            (self.handle, self.name, 'input_filename', handle)
        pystir.deleteDataHandle(handle)
    def set_zero_seg0_end_planes(self, flag):
        handle = pystir.charDataHandle(repr(flag))
        #pystir.cSTIR_setParameter\
        _setParameter\
            (self.handle, self.name, 'zero_seg0_end_planes', handle)
        pystir.deleteDataHandle(handle)
    def set_max_segment_num_to_process(self, n):
        set_int_par(self.handle, self.name, 'max_segment_num_to_process', n)
    def set_projector_pair_type(self, pp):
        #pystir.cSTIR_setParameter\
        _setParameter\
            (self.handle, self.name, 'projector_pair_type', pp.handle)

class Reconstruction:
    def set_output_filename_prefix(self, prefix):
        handle = pystir.charDataHandle(prefix)
        #pystir.cSTIR_setParameter\
        _setParameter\
            (self.handle, 'Reconstruction', 'output_filename_prefix', handle)
        pystir.deleteDataHandle(handle)

class IterativeReconstruction(Reconstruction):
    def __init__(self):
        self.handle = pystir.newDataHandle()
        self.name = 'IterativeReconstruction'
    def __del__(self):
        pystir.deleteDataHandle(self.handle)
    def set_num_subsets(self, n):
        set_int_par(self.handle, 'IterativeReconstruction', 'num_subsets', n)
    def get_num_subsets(self):
        return get_int_par(self.handle, 'IterativeReconstruction',\
                           'num_subsets')
    def set_start_subset_num(self, n):
        set_int_par(self.handle, 'IterativeReconstruction',\
                    'start_subset_num', n)
    def get_start_subset_num(self):
        return get_int_par(self.handle, 'IterativeReconstruction',\
                           'start_subset_num')
    def set_num_subiterations(self, n):
        set_int_par(self.handle, 'IterativeReconstruction',\
                    'num_subiterations', n)
    def get_num_subiterations(self):
        return get_int_par(self.handle, 'IterativeReconstruction',\
                           'num_subiterations')
    def set_start_subiteration_num(self, n):
        set_int_par(self.handle, 'IterativeReconstruction',\
                    'start_subiteration_num', n)
    def get_start_subiteration_num(self):
        return get_int_par(self.handle, 'IterativeReconstruction',\
                           'start_subiteration_num')
    def set_save_interval(self, n):
        set_int_par(self.handle, 'IterativeReconstruction',\
                    'save_interval', n)
    def set_inter_iteration_filter_interval(self, n):
        set_int_par(self.handle, 'IterativeReconstruction',\
                    'inter_iteration_filter_interval', n)
    def get_objective_function(self):
        obj_fun = GeneralisedObjectiveFunction()
        pystir.deleteDataHandle(obj_fun.handle)
        obj_fun.handle = pystir.cSTIR_parameter\
            (self.handle, 'IterativeReconstruction', 'objective_function')
        _check_status(obj_fun.handle)
        return obj_fun
    def set_objective_function(self, obj_fun):
        #pystir.cSTIR_setParameter\
        _setParameter\
            (self.handle, 'IterativeReconstruction',\
             'objective_function', obj_fun.handle)
    def set_inter_iteration_filter_type(self, filter_type):
        pystir.cSTIR_setParameter\
            (self.handle, 'IterativeReconstruction',\
             'inter_iteration_filter_type', filter_type.handle)
    def set_initial_estimate_filename(self, name):
        handle = pystir.charDataHandle(name)
        #pystir.cSTIR
        _setParameter\
            (self.handle, 'IterativeReconstruction', 'initial_estimate', handle)
        pystir.deleteDataHandle(handle)
    def set_up(self, image):
        handle = pystir.cSTIR_setupReconstruction(self.handle, image.handle)
        _check_status(handle)
        pystir.deleteDataHandle(handle)
    def reconstruct(self, image):
        handle = pystir.cSTIR_reconstruct(self.handle, image.handle)
        _check_status(handle)
        pystir.deleteDataHandle(handle)
    def set_subiteration_num(self, i):
        pystir.cSTIR_setSubiteration(self.handle, i)
    def subiteration_num(self):
        return pystir.cSTIR_getSubiteration(self.handle)
##        handle = pystir.cSTIR_getSubiteration(self.handle)
##        _check_status(handle)
##        iter = pystir.intDataFromHandle(handle)
##        pystir.deleteDataHandle(handle)
##        return iter;
    def update(self, image):
        handle = pystir.cSTIR_update(self.handle, image.handle)
        _check_status(handle)
        pystir.deleteDataHandle(handle)

class OSMAPOSLReconstruction3DFloat(IterativeReconstruction):
    def __init__(self, file = ''):
        self.name = 'OSMAPOSL'
        self.handle = pystir.cSTIR_newReconstruction(self.name, file)
        _check_status(self.handle)
    def __del__(self):
        pystir.deleteDataHandle(self.handle)
    def get_objective_function(self):
        obj_fun = PoissonLogLikelihoodWithLinearModelForMean3DFloat()
        pystir.deleteDataHandle(obj_fun.handle)
        obj_fun.handle = pystir.cSTIR_parameter\
            (self.handle, 'IterativeReconstruction', 'objective_function')
        _check_status(obj_fun.handle)
        return obj_fun
    def set_MAP_model(self, model):
        handle = pystir.charDataHandle(model)
        #pystir.cSTIR
        _setParameter\
            (self.handle, self.name, 'MAP_model', handle)
        pystir.deleteDataHandle(handle)

class OSSPSReconstruction3DFloat(IterativeReconstruction):
    def __init__(self, file = ''):
        self.name = 'OSSPS'
        self.handle = pystir.cSTIR_newReconstruction(self.name, file)
        _check_status(self.handle)
    def __del__(self):
        pystir.deleteDataHandle(self.handle)
    def set_relaxation_parameter(self, par):
        handle = pystir.floatDataHandle(par)
        #pystir.cSTIR
        _setParameter\
            (self.handle, self.name, 'relaxation_parameter', handle)
        pystir.deleteDataHandle(handle)
