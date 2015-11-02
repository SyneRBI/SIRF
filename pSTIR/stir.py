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
        errorMsg = \
            repr(msg) + ' exception thrown at line ' + \
            repr(line) + ' of ' + file
        raise error(errorMsg)

def _setParameter(hs, set, par, hv):
    h = pystir.cSTIR_setParameter(hs, set, par, hv)
    _check_status(h)
    pystir.deleteDataHandle(h)
def _set_int_par(handle, set, par, value):
    h = pystir.intDataHandle(value)
    _setParameter(handle, set, par, h)
    pystir.deleteDataHandle(h)
def _set_float_par(handle, set, par, value):
    h = pystir.floatDataHandle(value)
    _setParameter(handle, set, par, h)
    pystir.deleteDataHandle(h)
def _set_char_par(handle, set, par, value):
    h = pystir.charDataHandle(value)
    _setParameter(handle, set, par, h)
    pystir.deleteDataHandle(h)
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

class printerTo:
    def __init__(self, dest, channel = -1):
        if dest in {'stdout', 'stderr', 'cout', 'cerr'}:
            self.printer = pystir.newTextPrinter(dest)
            self.type = 0
        else:
            self.printer = pystir.newTextWriter(dest)
            self.type = 1
        pystir.openChannel(channel, self.printer)
        self.channel = channel
    def __del__(self):
        if self.type == 0:
            pystir.deleteTextPrinter(self.printer)
        else:
            pystir.deleteTextWriter(self.printer)
        pystir.closeChannel(self.channel)

class Shape:
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
        self.name = 'EllipsoidalCylinder'
        self.handle = pystir.cSTIR_newObject(self.name)
    def __del__(self):
        pystir.cSTIR_deleteObject(self.handle, 'Shape')
##        pystir.cSTIR_deleteObject(self.handle, self.name)
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
        self.ok = False
        self.name = 'Voxels'
        self.handle = pystir.cSTIR_voxels3DF\
                      (dim[0], dim[1], dim[2], \
                       vsize[0], vsize[1], vsize[2], \
                       origin[0], origin[1], origin[2])
        _check_status(self.handle)
        self.ok = True
    def __del__(self):
        if self.ok:
            pystir.cSTIR_deleteObject(self.handle, self.name)

class Image:
    def __init__(self, arg = None):
        self.empty = True
        self.voxels = None
        if isinstance(arg, str):
            self.handle = pystir.cSTIR_imageFromFile(arg)
            _check_status(self.handle)
            self.empty = False
        elif arg is not None:
            raise error\
                  ('wrong argument ' + repr(arg) + ' for Image constructor')
        self.name = 'Image'
        self.rimsize = -1
    def __del__(self):
        if not self.empty:
            pystir.cSTIR_deleteObject(self.handle, self.name)
        if self.voxels is not None:
            pystir.cSTIR_deleteObject(self.voxels, 'Voxels')
    def initialise(self, dim, vsize = (1, 1, 1), origin = (0, 0, 0)):
        if not self.empty:
            pystir.cSTIR_deleteObject(self.handle, self.name)
        if self.voxels is not None:
            pystir.cSTIR_deleteObject(self.voxels, 'Voxels')
        self.empty = True #temporary (blocks destructor's complains)
        self.voxels = pystir.cSTIR_voxels3DF\
                      (dim[0], dim[1], dim[2], \
                       vsize[0], vsize[1], vsize[2], \
                       origin[0], origin[1], origin[2])
        _check_status(self.voxels)
        self.handle = pystir.cSTIR_imageFromVoxels(self.voxels)
        _check_status(self.handle)
        self.empty = False
    def add_shape(self, shape, scale):
        handle = pystir.cSTIR_addShape\
                 (self.handle, self.voxels, shape.handle, scale)
        _check_status(handle)
        pystir.deleteDataHandle(handle)
    def read_from_file(self, filename):
        if not self.empty:
            pystir.cSTIR_deleteObject(self.handle, self.name)
        if self.voxels is not None:
            pystir.cSTIR_deleteObject(self.voxels, 'Voxels')
        self.empty = True
        self.handle = pystir.cSTIR_imageFromFile(filename)
        _check_status(self.handle)
        self.empty = False
    def diff_from(self, image):
        handle = pystir.cSTIR_imagesDifference\
                 (self.handle, image.handle, self.rimsize)
        _check_status(handle)
        diff = pystir.doubleDataFromHandle(handle)
        pystir.deleteDataHandle(handle)
        return diff
    def density(self):
        dim = numpy.ndarray((3,), dtype = numpy.int32)
        pystir.cSTIR_getImageDimensions(self.handle, dim.ctypes.data)
        nz = dim[0]
        ny = dim[1]
        nx = dim[2]
        density = numpy.ndarray((nz, ny, nx), dtype = numpy.float64)
        pystir.cSTIR_getImageData(self.handle, density.ctypes.data)
        return density

class TruncateToCylindricalFOVImageProcessor:
    def __init__(self):
        self.name = 'TruncateToCylindricalFOVImageProcessor'
        self.handle = pystir.cSTIR_newObject(self.name);
    def __del__(self):
        pystir.cSTIR_deleteObject(self.handle, 'DataProcessor')
##        pystir.cSTIR_deleteObject(self.handle, self.name)

class ProjectorsUsingMatrix:
    def __init__(self):
        self.name = 'ProjectorsUsingMatrix'
        self.handle = pystir.cSTIR_newObject(self.name);
        self.matrix = None
    def __del__(self):
        pystir.cSTIR_deleteObject(self.handle, 'Projectors')
##        pystir.cSTIR_deleteObject(self.handle, self.name)
    def set_matrix(self, matrix):
        _setParameter(self.handle, self.name, 'matrix_type', matrix.handle)
        self.matrix = matrix
    def get_matrix(self):
        if self.matrix is None:
            raise error(self.name + ': no matrix set')
        else:
            return self.matrix

class RayTracingMatrix:
    def __init__(self):
        self.name = 'RayTracingMatrix'
        self.handle = pystir.cSTIR_newObject(self.name);
        _check_status(self.handle)
    def __del__(self):
        pystir.cSTIR_deleteObject(self.handle, 'ProjMatrix')
##        pystir.cSTIR_deleteObject(self.handle, self.name)
    def set_num_tangential_LORs(self, value):
        _set_int_par(self.handle, self.name, 'num_tangential_LORs', value)
    def get_num_tangential_LORs(self):
        return _int_par(self.handle, self.name, 'num_tangential_LORs')

class GeneralisedPrior:
    def __init__(self):
        self.name_ = 'GeneralisedPrior'
    def set_penalisation_factor(self, value):
        _set_float_par(self.handle, self.name_, 'penalisation_factor', value)
    def get_penalisation_factor(self):
        return _float_par(self.handle, self.name_, 'penalisation_factor')
    def set_up(self):
        handle = pystir.cSTIR_setupObject('GeneralisedPrior', self.handle)
        _check_status(handle)
        pystir.deleteDataHandle(handle)

class QuadraticPrior(GeneralisedPrior):
    def __init__(self):
        GeneralisedPrior.__init__(self)
        self.name = 'QuadraticPrior'
        self.handle = pystir.cSTIR_newObject(self.name);
    def __del__(self):
        pystir.cSTIR_deleteObject(self.handle, 'Prior')
##        pystir.cSTIR_deleteObject(self.handle, self.name)

class GeneralisedObjectiveFunction:
    def __init__(self):
        self.prior = None
        self.name_ = 'GeneralisedObjectiveFunction'
    def set_prior(self, prior):
        _setParameter(self.handle, self.name_, 'prior', prior.handle)
        self.prior = prior
    def get_prior(self):
        if self.prior is None:
            raise error(self.name_, ': no prior set')
        else:
            return self.prior
    def set_up(self):
        handle = pystir.cSTIR_setupObject(self.name_, self.handle)
        _check_status(handle)
        pystir.deleteDataHandle(handle)

class PoissonLogLikelihoodWithLinearModelForMean\
      (GeneralisedObjectiveFunction):
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

class PoissonLogLikelihoodWithLinearModelForMeanAndProjData\
      (PoissonLogLikelihoodWithLinearModelForMean):
    def __init__(self):
        GeneralisedObjectiveFunction.__init__(self)
        self.name = 'PoissonLogLikelihoodWithLinearModelForMeanAndProjData'
        self.handle = pystir.cSTIR_newObject(self.name);
        self.projectors = None
    def __del__(self):
        pystir.cSTIR_deleteObject(self.handle, 'ObjectiveFunction')
##        pystir.cSTIR_deleteObject(self.handle, self.name)
    def set_input_filename(self, name):
        _set_char_par\
            (self.handle, self.name, 'input_filename', name)
    def set_zero_seg0_end_planes(self, flag):
        _set_char_par\
            (self.handle, self.name, 'zero_seg0_end_planes', repr(flag))
    def set_max_segment_num_to_process(self, n):
        _set_int_par(self.handle, self.name, 'max_segment_num_to_process', n)
    def set_projector_pair(self, pp):
        _setParameter\
            (self.handle, self.name, 'projector_pair_type', pp.handle)
        self.projectors = pp
    def get_projector_pair(self):
        if self.projectors is None:
            raise error(self.name + ': no projectors set')
        else:
            return self.projectors

class Reconstruction:
    def set_output_filename_prefix(self, prefix):
        _set_char_par\
            (self.handle, 'Reconstruction', 'output_filename_prefix', prefix)

class IterativeReconstruction(Reconstruction):
    def __init__(self):
        self.obj_fun = None
        self.filter = None
        self.name_ = 'IterativeReconstruction'
    def set_num_subsets(self, n):
        _set_int_par(self.handle, self.name_, 'num_subsets', n)
    def get_num_subsets(self):
        return _int_par(self.handle, self.name_, 'num_subsets')
    def set_start_subset_num(self, n):
        _set_int_par(self.handle, self.name_, 'start_subset_num', n)
    def get_start_subset_num(self):
        return _int_par(self.handle, self.name_, 'start_subset_num')
    def set_num_subiterations(self, n):
        _set_int_par(self.handle, self.name_, 'num_subiterations', n)
    def get_num_subiterations(self):
        return _int_par(self.handle, self.name_, 'num_subiterations')
    def set_start_subiteration_num(self, n):
        _set_int_par(self.handle, self.name_, 'start_subiteration_num', n)
    def get_start_subiteration_num(self):
        return _int_par(self.handle, self.name_, 'start_subiteration_num')
    def set_subiteration_num(self, iter):
        _set_int_par(self.handle, self.name_, 'subiteration_num', iter)
    def get_subiteration_num(self):
        return _int_par(self.handle, self.name_, 'subiteration_num')
    def set_save_interval(self, n):
        _set_int_par(self.handle, self.name_, 'save_interval', n)
    def set_inter_iteration_filter_interval(self, n):
        _set_int_par(self.handle, self.name_, 'inter_iteration_filter_interval', n)
    def set_objective_function(self, obj):
        _setParameter(self.handle, self.name_, 'objective_function', obj.handle)
        self.obj_fun = obj
    def get_objective_function(self):
        if self.obj_fun is None:
            raise error(self.name_ + ': no objective function set')
        else:
            return self.obj_fun
    def set_inter_iteration_filter(self, f):
        pystir.cSTIR_setParameter\
            (self.handle, self.name_, 'inter_iteration_filter_type', f.handle)
    def get_inter_iteration_filter(self):
        if self.filter is None:
            raise error(self.name_ + ': no filter set')
        else:
            return self.filter
    def set_up(self, image):
        handle = pystir.cSTIR_setupReconstruction(self.handle, image.handle)
        _check_status(handle)
        pystir.deleteDataHandle(handle)
    def reconstruct(self, image):
        handle = pystir.cSTIR_reconstruct(self.handle, image.handle)
        _check_status(handle)
        pystir.deleteDataHandle(handle)
    def update(self, image):
        handle = pystir.cSTIR_update(self.handle, image.handle)
        _check_status(handle)
        pystir.deleteDataHandle(handle)

class OSMAPOSLReconstruction(IterativeReconstruction):
    def __init__(self, file = ''):
        super(OSMAPOSLReconstruction, self).__init__()
        self.ok = False
        self.name = 'OSMAPOSL'
        self.handle = pystir.cSTIR_newReconstruction(self.name, file)
        _check_status(self.handle)
        self.ok = True
    def __del__(self):
        if self.ok:
            pystir.cSTIR_deleteObject(self.handle, 'Reconstruction')
##            pystir.cSTIR_deleteReconstruction(self.handle)
    def set_MAP_model(self, model):
        _set_char_par\
            (self.handle, self.name, 'MAP_model', model)

class OSSPSReconstruction(IterativeReconstruction):
    def __init__(self, file = ''):
        super(OSSPSReconstruction, self).__init__()
        self.ok = False
        self.name = 'OSSPS'
        self.handle = pystir.cSTIR_newReconstruction(self.name, file)
        _check_status(self.handle)
        self.ok = True
    def __del__(self):
        if self.ok:
            pystir.cSTIR_deleteObject(self.handle, 'Reconstruction')
##            pystir.cSTIR_deleteReconstruction(self.handle)
    def set_relaxation_parameter(self, value):
        _set_float_par\
            (self.handle, self.name, 'relaxation_parameter', value)
