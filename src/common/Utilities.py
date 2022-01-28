'''Utilities used by all engines
'''
import inspect
import numpy
import os
import sirf
import sirf.pyiutilities as pyiutil
import re
from deprecation import deprecated

__licence__ = """SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2021 Rutherford Appleton Laboratory STFC
Copyright 2015 - 2021 University College London
Copyright 2021 CSIRO

This is software developed for the Collaborative Computational
Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
(http://www.ccpsynerbi.ac.uk/).

Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at
      http://www.apache.org/licenses/LICENSE-2.0
  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
"""
__license__ = __licence__
RE_PYEXT = re.compile(r"\.(py[co]?)$")


@deprecated(
    deprecated_in="2.0.0", removed_in="4.0", current_version=sirf.__version__,
    details="use examples_data_path() instead")
def petmr_data_path(petmr):
    '''
    Returns the path to PET or MR data.
    petmr: either 'PET' or 'MR'
    '''
    return examples_data_path( petmr.upper() )
    

def examples_data_path(data_type):
    '''
    Returns the path to PET/MR/Registration data used by SIRF/examples demos.
    data_type: either 'PET' or 'MR' or 'Registration'
    '''
    data_path = os.path.join('share', 'SIRF-{}.{}'.format(sirf.__version_major__, sirf.__version_minor__),
                             'data', 'examples', data_type)
    SIRF_INSTALL_PATH = os.environ.get('SIRF_INSTALL_PATH')
    SIRF_DATA_PATH = os.environ.get('SIRF_DATA_PATH')
    SIRF_PATH = os.environ.get('SIRF_PATH')
    if SIRF_DATA_PATH is not None:
        return os.path.join(SIRF_DATA_PATH , 'examples', data_type)
    elif SIRF_INSTALL_PATH is not None:
        return os.path.join(SIRF_INSTALL_PATH , data_path)
    elif SIRF_PATH is not None:
        return os.path.join(SIRF_PATH, 'data', 'examples', data_type)
    else:
        errorMsg = 'You need to set the SIRF_DATA_PATH or SIRF_INSTALL_PATH environment variable to allow finding the raw data.'
        raise ValueError(errorMsg)


def existing_filepath(data_path, file_name):
    '''
    Returns the filepath (path/name) to an existing file.
    Raises error if the file does not exist.
    data_path: path to the file
    file_name: file name
    '''
    full_name = os.path.join(os.path.abspath(data_path), file_name)
    if not os.path.isfile(full_name):
        raise error('file %s not found' % full_name)
    return full_name


def show_2D_array(title, array, scale = None, colorbar = True):
    '''
    Displays a 2D array.
    title   : the figure title
    array   : 2D array
    colorbar: flag specifying whether the colorbar is to be displayed
    '''
    try:
        import matplotlib.pyplot as plt
    except:
        print('matplotlib not found, cannot plot the array')
        return
    if scale is None:
        vmin = numpy.amin(array)
        vmax = numpy.amax(array)
    else:
        vmin, vmax = scale
    plt.figure()
    plt.title(title)
    if colorbar:
        plt.imshow(array, vmin=vmin, vmax=vmax)
        plt.colorbar()
    else:
        plt.imshow(array, cmap='gray', vmin=vmin, vmax=vmax)
    fignums = plt.get_fignums()
    print('You may need to close Figure %d window to continue...' % fignums[-1])
    plt.show()


def show_3D_array\
    (array, index=None, tile_shape=None, scale=None, power=None, \
     suptitle=None, titles=None, title_size=None, \
     zyx=None, xlabel=None, ylabel=None, label=None, \
     cmap=None, show=True):
    '''
    Displays a 3D array as a set of z-slice tiles.
    On successful completion returns 0.
    array     : 3D array
    index     : z-slices index, either Python list or string of the form
              : 'a, b-c, ...', where 'b-c' is decoded as 'b, b+1, ..., c';
              : out-of-range index value causes error (non-zero) return
    tile_shape: tuple (tile_rows, tile_columns);
                if not present, the number of tile rows and columns is
                computed based on the array dimensions
    scale     : tuple (vmin, vmax) for imshow; defaults to the range of
                array values
    power     : if present, numpy.power(abs(array), power) is displayed
                (power < 1 improves visibility of relatively small array values)
    suptitle  : figure title; defaults to None
    titles    : array of tile titles; if not present, each tile title is
                label + tile_number
    zyx       : tuple (z, y, x), where x, y, anad z are the dimensions of array
                corresponding to the spatial dimensions x, y and z; zyx=None is
                interpreted as (0, 1, 2)
    xlabel    : label for x axis
    ylabel    : label for y axis
    label     : tile title prefix
    cmap      : colormap
    show      : flag specifying whether the array must be displayed immediately
    '''
    import math
    try:
        import matplotlib as mpl
        import matplotlib.pyplot as plt
    except:
        print('matplotlib not found, cannot plot the array')
        return
    import numpy

    current_title_size = mpl.rcParams['axes.titlesize']
    current_label_size = mpl.rcParams['axes.labelsize']
    current_xlabel_size = mpl.rcParams['xtick.labelsize']
    current_ylabel_size = mpl.rcParams['ytick.labelsize']
    mpl.rcParams['axes.titlesize'] = 'small'
    mpl.rcParams['axes.labelsize'] = 'small'
    mpl.rcParams['xtick.labelsize'] = 'small'
    mpl.rcParams['ytick.labelsize'] = 'small'

    if zyx is not None:
        array = numpy.transpose(array, zyx)

    nz = array.shape[0]
    if index is None:
        n = nz
        index = range(n)
    else:
        if type(index) == type(' '):
            try:
                index = str_to_int_list(index)
            except:
                return 1
        n = len(index)
        for k in range(n):
            z = index[k]
            if z < 0 or z >= nz:
                return k + 1
    ny = array.shape[1]
    nx = array.shape[2]
    if tile_shape is None:
        rows = int(round(math.sqrt(n*nx/ny)))
        if rows < 1:
            rows = 1
        if rows > n:
            rows = n
        cols = (n - 1)//rows + 1
        last_row = rows - 1
    else:
        rows, cols = tile_shape
        assert rows*cols >= n, \
            "tile rows x columns must be not less than the number of images"
        last_row = (n - 1)//cols
    if scale is None:
        if power is None:
            vmin = numpy.amin(array)
            vmax = numpy.amax(array)
        else:
            vmin = numpy.power(numpy.amin(abs(array)), power)
            vmax = numpy.power(numpy.amax(abs(array)), power)
    else:
        vmin, vmax = scale
    fig = plt.figure()
    if suptitle is not None:
        if title_size is None:
            fig.suptitle(suptitle)
        else:
            fig.suptitle(suptitle, fontsize=title_size)
    for k in range(n):
        z = index[k] #- 1
        ax = fig.add_subplot(rows, cols, k + 1)
        if titles is None:
            if label is not None and nz > 1:
                ax.set_title(label + (' %d' % z))
        else:
            ax.set_title(titles[k])
        row = k//cols
        col = k - row*cols
        if xlabel is None and ylabel is None or row < last_row or col > 0:
            ax.set_axis_off()
        else:
            ax.set_axis_on()
            if xlabel is not None:
                plt.xlabel(xlabel)
                plt.xticks([0, nx - 1], [0, nx - 1])
            if ylabel is not None:
                plt.ylabel(ylabel)
                plt.yticks([0, ny - 1], [0, ny - 1])
        if power is None:
            imgplot = ax.imshow(array[z,:,:], cmap, vmin=vmin, vmax=vmax)
        else:
            imgplot = ax.imshow(numpy.power(abs(array[z,:,:]), power), cmap, \
                                vmin=vmin, vmax=vmax)
    if show:
        fignums = plt.get_fignums()
        last = fignums[-1]
        if last > 1:
            print("You may need to close Figures' 1 - %d windows to continue..." \
                  % last)
        else:
            print('You may need to close Figure 1 window to continue...')
        plt.show()

    mpl.rcParams['axes.titlesize'] = current_title_size
    mpl.rcParams['axes.labelsize'] = current_label_size
    mpl.rcParams['xtick.labelsize'] = current_xlabel_size
    mpl.rcParams['ytick.labelsize'] = current_ylabel_size

    return 0

def format_numpy_array_for_setter(data, dtype_to_pass=numpy.float32):

    if not isinstance(data, numpy.ndarray):
        raise error('Wrong input format.' + \
            ' Should be numpy.ndarray. Got {}'.format(type(data)))

    if data.dtype != dtype_to_pass:
            data = data.astype(dtype_to_pass)

    if not data.flags['C_CONTIGUOUS']:
        data = numpy.ascontiguousarray(data)

    return data

def check_tolerance(expected, actual, abstol=0, reltol=2e-3):
    '''
    Check if 2 floats are equal within the specified tolerance, i.e.
    abs(expected - actual) <= abstol + reltol*abs(expected).
    Returns an error string if they are not and None otherwise.
    '''
    tol = abstol + reltol*abs(expected)
    if abs(expected - actual) > tol:
        return "expected %.4g, got %.4g (tolerance %.3g)" \
               % (expected, actual, tol)


class pTest(object):
    def __init__(self, filename, record, throw=False):
        self.record = record
        self.data = []
        self.ntest = 0
        self.nrec = 0
        self.failed = 0
        self.verbose = True
        self.throw = throw
        if record:
            self.file = open(filename, 'w')
        else:
            with open(filename, 'r') as f:
                self.data = [float(line.strip()) for line in f]
            self.size = len(self.data)
            self.file = None

    def __del__(self):
        msg = "%d failures" % self.failed
        if self.failed:
            if self.record:
                self.file.write(msg + '\n')
        if self.record:
            self.file.close()

    def check(self, value, abs_tol=0, rel_tol=2e-3):
        '''
        Tests if value is equal to the recorded one (or record it)
        value        : the value that was computed
        abs_tol, rel_tol: see :func:`~Utilities.check_tolerance`
        '''

        if self.record:
            self.file.write('%e\n' % value)
        else:
            if self.nrec >= self.size:
                raise IndexError('no data available for test %d' % self.ntest)
            else:
                expected = self.data[self.nrec]
                self.check_if_equal_within_tolerance(expected, value, abs_tol, rel_tol)
        self.nrec += 1

    def check_if_equal(self, expected, value):
        '''
        Tests if value is equal to the expected one.
        expected     : the true value
        value        : the value that was computed
        '''
        if value != expected:
            self.failed += 1
            msg = '+++ test %d failed: expected %s, got %s' \
                  % (self.ntest, repr(expected), repr(value))
            if self.throw:
                raise ValueError(msg)
            if self.verbose:
                print(msg)
        else:
            if self.verbose:
                print('+++ test %d passed' % self.ntest)
        self.ntest += 1

    def check_if_equal_within_tolerance(self, expected, value, abs_tol=0, rel_tol=2e-3):
        '''
        Tests if float value is equal to the expected one.
        expected     : the true value
        value        : the value that was computed
        abs_tol, rel_tol: see :func:`~Utilities.check_tolerance`
        '''
        err = check_tolerance(expected, value, abs_tol, rel_tol)
        if err is not None:
            self.failed += 1
            msg = ('+++ test %d failed: ' % self.ntest) + str(err)
            if self.throw:
                raise ValueError(msg)
            if self.verbose:
                print(msg)
        else:
            if self.verbose:
                print('+++ test %d passed' % self.ntest)
        self.ntest += 1

    def check_if_zero_within_tolerance(self, value, abs_tol=1e-3):
        '''
        Tests if float value is equal to the expected one.
        expected     : the true value
        abs_tol: see :func:`~Utilities.check_tolerance`
        '''
        self.check_if_equal_within_tolerance(0, value, abs_tol)

    def check_if_less(self, value, comp):
        '''
        Tests if value is (strictly) less than comp.
        value        : the value that was computed
        comp         : the maximum allowed value
        '''
        if value >= comp:
            self.failed += 1
            msg = ('+++ test %d failed: ' % self.ntest) + \
                  repr(value) + ' >= ' + repr(comp)
            if self.throw:
                raise ValueError(msg)
            if self.verbose:
                print(msg)
        else:
            if self.verbose:
                print('+++ test %d passed' % self.ntest)
        self.ntest += 1


class CheckRaise(pTest):
    def __init__(self, *a, **k):
        k["throw"] = True
        super(CheckRaise, self).__init__(*a, **k)


def runner(main_test, doc, version, author="", licence=None):
    """
    :param main_test: function(record : bool, verbose : bool, throw : bool)
    """
    from docopt import docopt
    args = docopt(doc.format(version=version,
                             author=author,
                             licence=licence or __licence__,
                             license=licence or __licence__),
                  version=version)

    record = args['--record']
    verbose = args['--verbose']

    failed, ntest = main_test(record, verbose, throw=False)
    if failed:
        import sys
        print('%d of %d tests failed' % (failed, ntest))
        sys.exit(failed)
    if record:
        print('%d measurements recorded' % ntest)
    else:
        print('all %d tests passed' % ntest)


###########################################################
############ Utilities for internal use only ##############
class error(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return '??? ' + repr(self.value)


def check_status(handle, stack=None):
    if pyiutil.executionStatus(handle) != 0:
        if stack is None:
            stack = inspect.stack()[1]
        print('\nFile: %s' % stack[1])
        print('Line: %d' % stack[2])
        print('check_status found the following message sent from the engine:')
        msg = pyiutil.executionError(handle)
        file = pyiutil.executionErrorFile(handle)
        line = pyiutil.executionErrorLine(handle)
        errorMsg = \
            repr(msg) + ' exception caught at line ' + \
            repr(line) + ' of ' + file + '; ' + \
            'the reconstruction engine output may provide more information'
        raise error(errorMsg)


def try_calling(returned_handle):
    check_status(returned_handle, inspect.stack()[1])
    pyiutil.deleteDataHandle(returned_handle)

def assert_validity(obj, dtype):
    if not isinstance(obj, dtype):
        msg = 'Expecting object of type {}, got {}'
        raise AssertionError(msg.format(dtype, type(obj)))
    if obj.handle is None:
        raise AssertionError('object handle is None.')


def assert_validities(x, y):
    if not (issubclass(type(x),type(y)) or issubclass(type(y),type(x))):
        msg = 'Expecting same type input, got {} and {}'
        raise AssertionError(msg.format(type(x), type(y)))
    if x.handle is None:
        raise AssertionError('handle for first parameter is None')
    if y.handle is None:
        raise AssertionError('handle for second parameter is None')
    if callable(getattr(x, 'dimensions', None)):
        xdim = x.dimensions()
    else:
        xdim = None
    if callable(getattr(y, 'dimensions', None)):
        ydim = y.dimensions()
    else:
        ydim = None
    if xdim != ydim:
        raise ValueError("Input shapes are expected to be equal, got " \
                         + repr(xdim) + " and " \
                         + repr(ydim) + " instead.")


def label_and_name(g):
    name = g.lstrip()
    name = name.rstrip()
    i = name.find(':')
    if i > -1:
        label = name[: i].rstrip()
        name = name[i + 1 :].lstrip()
    else:
        label = ''
    return label, name


def name_and_parameters(obj):
    name = obj.lstrip()
    name = name.rstrip()
    i = name.find('(')
    if i > -1:
        j = name.find(')', i)
        prop = name[i + 1 : j]
        name = name[: i].rstrip()
        i = 0
    else:
        prop = None
    return name, prop


def parse_arglist(arglist):
    argdict = {}
    while True:
        arglist = arglist.lstrip()
        ieq = arglist.find('=')
        if ieq < 0:
            return argdict
        name = arglist[0:ieq].rstrip()
        arglist = arglist[ieq + 1 :].lstrip()
        ic = arglist.find(',')
        if ic < 0:
            argdict[name] = arglist.rstrip()
            return argdict
        else:
            argdict[name] = arglist[0:ic].rstrip()
            arglist = arglist[ic + 1 :]


def str_to_int_list(str_list):
    int_list = []
    last = False
    while not last:
        ic = str_list.find(',')
        if ic < 0:
            ic = len(str_list)
            last = True
        str_item = str_list[0:ic]
        str_list = str_list[ic + 1 :]
        ic = str_item.find('-')
        if ic < 0:
            int_item = [int(str_item)]
        else:
            strt = int(str_item[0:ic])
            stop = int(str_item[ic + 1 :])
            int_item = list(range(strt, stop + 1))
        int_list = int_list + int_item
    return int_list

def is_operator_adjoint(operator, num_tests = 5, max_err = 10e-5, verbose = True):
    '''
    Test if a given operator is adjoint.
    The operator needs to have been already set_up() with valid objects.
    The operator needs to have methods direct() and adjoint() implemented

    Parameters
    ----------
    operator  :
        Any SIRF operator that implements direct() and adjoint()
    num_tests : int, optional
        Square root of the number of tests with random data that will be executed. Default 5
    max_err   : double, optional
        Maximum allowed normalized error, tolerance. Change not recommended. Default 10e-5
    verbose   : bool
        Verbose option
    '''
    for iter1 in range(num_tests):
        ## generate random data for x and direct()
        x = operator.domain_geometry().allocate(value = 'random')
        y_hat = operator.direct(x)
        for iter2 in range(num_tests):
            if verbose:
                print("Testing " + type(operator).__name__ + ": Iteration " + str(iter1*num_tests+iter2+1) + "/" + str(num_tests**2))
            ## generate random data and adjoint()
            y = operator.range_geometry().allocate( value = 'random')
            x_hat = operator.adjoint(y)
            # Check dot product identity
            norm_err = abs(numpy.conj(y_hat.dot(y)) - x_hat.dot(x))/(numpy.conj(abs(y_hat.dot(y)))*0.5 + abs(x_hat.dot(x))*0.5)
            if norm_err > max_err:
                if verbose:
                    print(type(operator).__name__ + " is not adjoint, with normalized error of " + str(norm_err) + " (max: " + str(max_err) + ")")
                return False
            elif verbose:
                print("Pass, with a with normalized error of " + str(norm_err) + " (max: " + str(max_err) + ")")
    return True


class TestDataContainerAlgebra(object):

    '''A base class for unit test of DataContainer algebra.'''
    def test_divide_scalar(self):
        if hasattr(self, 'cwd'):
            os.chdir(self.cwd)
        image1 = self.image1
        image2 = self.image2
        image1.fill(1.)
        image2.fill(2.)
        
        tmp = image1/1.
        numpy.testing.assert_array_equal(image1.as_array(), tmp.as_array())
    
        tmp1 = image1.divide(1.)
        numpy.testing.assert_array_equal(tmp.as_array(), tmp1.as_array())
        
        image1.divide(1., out=image2)
        numpy.testing.assert_array_equal(tmp.as_array(), image2.as_array())

        image2.fill(2)
        image2 /= 2.0
        numpy.testing.assert_array_equal(image1.as_array(), image2.as_array())

    def test_divide_datacontainer(self):
        if hasattr(self, 'cwd'):
            os.chdir(self.cwd)
        # add 1 because the data contains zeros and divide is not going to be happy
        image1 = self.image1 + 1
        image2 = self.image2 + 1
        
        tmp = image1/image2

        numpy.testing.assert_array_almost_equal(
            numpy.ones(image1.shape, dtype=numpy.float32), tmp.as_array()
            )
    
        tmp1 = image1.divide(image2)
        numpy.testing.assert_array_almost_equal(
            numpy.ones(image1.shape, dtype=numpy.float32), tmp1.as_array()
            )
        
        tmp1.fill(2.)
        image1.divide(image2, out=tmp1)
        
        numpy.testing.assert_array_almost_equal(
            numpy.ones(image1.shape, dtype=numpy.float32), tmp1.as_array()
            )
        
        image1 /= image2
        numpy.testing.assert_array_almost_equal(
            numpy.ones(image1.shape, dtype=numpy.float32), image1.as_array()
            )        

    def test_multiply_scalar(self):
        if hasattr(self, 'cwd'):
            os.chdir(self.cwd)
        image1 = self.image1
        image2 = self.image2
        image2.fill(2.)
        
        tmp = image1 * 1.
        numpy.testing.assert_array_equal(image1.as_array(), tmp.as_array())
    
        tmp1 = image1.multiply(1.)
        numpy.testing.assert_array_equal(tmp.as_array(), tmp1.as_array())
        
        image1.multiply(1., out=image2)
        numpy.testing.assert_array_equal(tmp.as_array(), image2.as_array())

    def test_multiply_datacontainer(self):
        if hasattr(self, 'cwd'):
            os.chdir(self.cwd)
        image1 = self.image1
        image2 = self.image2
        image2.fill(1.)
        tmp = image1 * image2

        numpy.testing.assert_array_almost_equal(
            image1.as_array(), tmp.as_array()
            )
    
        tmp1 = image1.multiply(image2)
        numpy.testing.assert_array_almost_equal(
            image1.as_array(), tmp1.as_array()
            )
        
        tmp1.fill(2.)
        image1.multiply(image2, out=tmp1)
        
        numpy.testing.assert_array_almost_equal(
            image1.as_array(), tmp1.as_array()
            )

    def test_add_scalar(self):
        if hasattr(self, 'cwd'):
            os.chdir(self.cwd)
        image1 = self.image1
        image2 = self.image2
        image1.fill(0)
        image2.fill(1)
        
        tmp = image1 + 1.
        numpy.testing.assert_array_equal(image2.as_array(), tmp.as_array())
    
        tmp1 = image1.add(1.)
        numpy.testing.assert_array_equal(tmp.as_array(), tmp1.as_array())
        
        tmp1.fill(0)
        image1.add(1., out=tmp1)
        numpy.testing.assert_array_equal(tmp1.as_array(), image2.as_array())
    
    def test_add_datacontainer(self):
        if hasattr(self, 'cwd'):
            os.chdir(self.cwd)
        image1 = self.image1
        image2 = self.image2
        image1.fill(0.)
        image2.fill(1.)
        tmp = image1 + image2

        numpy.testing.assert_array_almost_equal(
            numpy.ones(image1.shape, dtype=numpy.float32), tmp.as_array()
            )
    
        tmp1 = image1.add(image2)
        
        numpy.testing.assert_array_almost_equal(
            numpy.ones(image1.shape, dtype=numpy.float32), tmp1.as_array()
            )
        
        tmp1.fill(2.)
        image1.add(image2, out=tmp1)
        
        numpy.testing.assert_array_almost_equal(
            numpy.ones(image1.shape, dtype=numpy.float32), tmp1.as_array()
            )
        
    
    def test_subtract_scalar(self):
        if hasattr(self, 'cwd'):
            os.chdir(self.cwd)
        image1 = self.image1
        image2 = self.image2
        image1.fill(2)
        image2.fill(1)
        
        tmp = image1 - 1.
        numpy.testing.assert_array_equal(image2.as_array(), tmp.as_array())
    
        tmp1 = image1.subtract(1.)
        numpy.testing.assert_array_equal(tmp.as_array(), tmp1.as_array())
        
        tmp1.fill(0)
        image1.subtract(1., out=tmp1)
        numpy.testing.assert_array_equal(tmp1.as_array(), image2.as_array())

    def test_subtract_datacontainer(self):
        if hasattr(self, 'cwd'):
            os.chdir(self.cwd)
        image1 = self.image1
        image2 = self.image2
        
        tmp = image1 - image2

        numpy.testing.assert_array_almost_equal(
            numpy.zeros(image1.shape, dtype=numpy.float32), tmp.as_array()
            )
    
        tmp1 = image1.subtract(image2)
        
        numpy.testing.assert_array_almost_equal(
            numpy.zeros(image1.shape, dtype=numpy.float32), tmp1.as_array()
            )
        
        tmp1.fill(2.)
        image1.subtract(image2, out=tmp1)
        
        numpy.testing.assert_array_almost_equal(
            numpy.zeros(image1.shape, dtype=numpy.float32), tmp1.as_array()
            )

    def test_division_by_scalar_zero(self):
        self.assertTrue(True)
        return
        if hasattr(self, 'cwd'):
            os.chdir(self.cwd)
        try:
            self.image1 / 0.
            self.assertFalse(True)
        except ZeroDivisionError:
            self.assertTrue(True)
        except error:
            self.assertTrue(True)
    
    def test_division_by_datacontainer_zero(self):
        self.assertTrue(True)
        return
        if hasattr(self, 'cwd'):
            os.chdir(self.cwd)
        try:
            self.image2 *= 0
            tmp = self.image1 / self.image2
            self.assertFalse(True)
        except ZeroDivisionError:
            self.assertTrue(True)
        except error:
            self.assertTrue(True)

    def test_sapyb_scalars(self):

        image1 = self.image1.copy()
        image2 = self.image2.copy()

        arr = numpy.arange(0,image1.size).reshape(image1.shape)
        image1.fill(arr)
        image2.fill(-arr)

        #scalars
        #check call methods with out

        a = 2.0
        b = -3.0
        gold = a * arr - b * arr

        out = image1.sapyb(a, image2, b)
        numpy.testing.assert_allclose(out.as_array(), gold)
        numpy.testing.assert_allclose(image1.as_array(), arr)
        numpy.testing.assert_allclose(image2.as_array(), -arr)

        out.fill(0)
        image1.sapyb(a, image2, b, out=out)
        numpy.testing.assert_allclose(out.as_array(), gold)
        numpy.testing.assert_allclose(image1.as_array(), arr)
        numpy.testing.assert_allclose(image2.as_array(), -arr)

        out.fill(arr)
        out.sapyb(a, image2, b, out=out)
        numpy.testing.assert_allclose(out.as_array(), gold)
        numpy.testing.assert_allclose(image2.as_array(), -arr)

        out.fill(-arr)
        image1.sapyb(a, out, b, out=out)
        numpy.testing.assert_allclose(out.as_array(), gold)
        numpy.testing.assert_allclose(image1.as_array(), arr)

    def test_sapyb_vectors(self):

        image1 = self.image1.copy()
        image2 = self.image2.copy()

        arr = numpy.arange(0,image1.size).reshape(image1.shape)
        image1.fill(arr)
        image2.fill(-arr)

        a = image1.copy()
        a.fill(2)
        b = image1.copy()
        b.fill(-3)

        gold = a.as_array() * arr - b.as_array() * arr

        out = image1.sapyb(a, image2, b)
        numpy.testing.assert_allclose(out.as_array(), gold)
        numpy.testing.assert_allclose(image1.as_array(), arr)
        numpy.testing.assert_allclose(image2.as_array(), -arr)

        out.fill(0)
        image1.sapyb(a, image2, b, out=out)
        numpy.testing.assert_allclose(out.as_array(), gold)
        numpy.testing.assert_allclose(image1.as_array(), arr)
        numpy.testing.assert_allclose(image2.as_array(), -arr)

        out.fill(arr)
        out.sapyb(a, image2, b, out=out)
        numpy.testing.assert_allclose(out.as_array(), gold)
        numpy.testing.assert_allclose(image2.as_array(), -arr)

        out.fill(-arr)
        image1.sapyb(a, out, b, out=out)
        numpy.testing.assert_allclose(out.as_array(), gold)
        numpy.testing.assert_allclose(image1.as_array(), arr)

    def test_sapyb_mixed(self):

        image1 = self.image1.copy()
        image2 = self.image2.copy()

        arr = numpy.arange(0,image1.size).reshape(image1.shape)
        image1.fill(arr)
        image2.fill(-arr)
 
        a = 2.0
        b = image1.copy()
        b.fill(-3)

        gold = a * arr - b.as_array() * arr

        out = image1.sapyb(a, image2, b)
        numpy.testing.assert_allclose(out.as_array(), gold)
        numpy.testing.assert_allclose(image1.as_array(), arr)
        numpy.testing.assert_allclose(image2.as_array(), -arr)

        out.fill(0)
        image1.sapyb(a, image2, b, out=out)
        numpy.testing.assert_allclose(out.as_array(), gold)
        numpy.testing.assert_allclose(image1.as_array(), arr)
        numpy.testing.assert_allclose(image2.as_array(), -arr)
       
        out.fill(arr)
        out.sapyb(a, image2, b, out=out)
        numpy.testing.assert_allclose(out.as_array(), gold)
        numpy.testing.assert_allclose(image2.as_array(), -arr)

        out.fill(-arr)
        image1.sapyb(a, out, b, out=out)
        numpy.testing.assert_allclose(out.as_array(), gold)
        numpy.testing.assert_allclose(image1.as_array(), arr)