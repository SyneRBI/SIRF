'''Utilities used by all engines
'''
import inspect
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy
import os
import pyiutilities as pyiutil
import re

__licence__ = """CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC

This is software developed for the Collaborative Computational
Project in Positron Emission Tomography and Magnetic Resonance imaging
(http://www.ccppetmr.ac.uk/).

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


def petmr_data_path(petmr):
    '''
    Returns the path to PET or MR data.
    petmr: either 'pet' or 'mr' (case-insensitive)
    '''
    data_path = '/data/examples/' + petmr.upper()
    SIRF_PATH = os.environ.get('SIRF_PATH')
    if SIRF_PATH is not None:
        return SIRF_PATH + data_path
    else:
        errorMsg = 'You need to set the SIRF_PATH environment variable to allow finding the raw data.'
        raise error(errorMsg)


def existing_filepath(data_path, file_name):
    '''
    Returns the filepath (path/name) to an existing file.
    Raises error if the file does not exist.
    data_path: path to the file
    file_name: file name
    '''
    full_name = data_path + '/' + file_name
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
    if scale is None:
        vmin = numpy.amin(array)
        vmax = numpy.amax(array)
    else:
        vmin, vmax = scale
    plt.figure()
    plt.title(title)
    plt.imshow(array, vmin = vmin, vmax = vmax)
    if colorbar:
        plt.colorbar()
    fignums = plt.get_fignums()
    print('Close Figure %d window to continue...' % fignums[-1])
    plt.show()


def show_3D_array\
    (array, index = None, tile_shape = None, scale = None, power = None, \
     suptitle = None, titles = None, \
     xlabel = None, ylabel = None, label = None, \
     title_size = None, \
     show = True):
    '''
    Displays a 3D array as a set of z-slice tiles.
    On successful completion returns 0.
    array     : 3D array
    index     : z-slices index (1-based), either Python list or string of the form
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
    xlabel    : label for x axis
    ylabel    : label for y axis
    label     : tile title prefix
    show      : flag specifying whether the array must be displayed immediately
    '''
    import math
    import numpy
    current_title_size = mpl.rcParams['axes.titlesize']
    current_label_size = mpl.rcParams['axes.labelsize']
    current_xlabel_size = mpl.rcParams['xtick.labelsize']
    current_ylabel_size = mpl.rcParams['ytick.labelsize']
    mpl.rcParams['axes.titlesize'] = 'small'
    mpl.rcParams['axes.labelsize'] = 'small'
    mpl.rcParams['xtick.labelsize'] = 'small'
    mpl.rcParams['ytick.labelsize'] = 'small'
    nz = array.shape[0]
    if index is None:
        n = nz
        index = range(1, n + 1)
    else:
        if type(index) == type(' '):
            try:
                index = str_to_int_list(index)
            except:
                print('incorrect input')
                return 0
        n = len(index)
        for k in range(n):
            z = index[k]
            if z < 1 or z > nz:
                return k + 1
    if tile_shape is None:
        ny = array.shape[1]
        nx = array.shape[2]
        rows = int(round(math.sqrt(n*nx/ny)))
        if rows < 1:
            rows = 1
        if rows > n:
            rows = n
        cols = (n - 1)//rows + 1
    else:
        rows, cols = tile_shape
        assert rows*cols >= array.shape[0],\
                "tile rows x columns must equal the 3rd dim extent of array"
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
            fig.suptitle(suptitle, fontsize = title_size)
    for k in range(n):
        z = index[k] - 1
        ax = fig.add_subplot(rows, cols, k + 1)
        if titles is None:
            if label is not None and nz > 1:
                ax.set_title(label + (' %d' % (z + 1)))
        else:
            ax.set_title(titles[k])
        row = k//cols
        col = k - row*cols
        if xlabel is None and ylabel is None or row < rows - 1 or col > 0:
            ax.set_axis_off()
        else:
            ax.set_axis_on()
            if xlabel is not None:
                plt.xlabel(xlabel)
                plt.xticks([0, nx - 1], [1, nx])
            if ylabel is not None:
                plt.ylabel(ylabel)
                plt.yticks([0, ny - 1], [1, ny])
        if power is None:
            imgplot = ax.imshow(array[z,:,:], vmin = vmin, vmax = vmax)
        else:
            imgplot = ax.imshow(numpy.power(abs(array[z,:,:]), power), \
                                vmin = vmin, vmax = vmax)
    if show:
        fignums = plt.get_fignums()
        last = fignums[-1]
        if last > 1:
            print("Close Figures' 1 - %d windows to continue..." % last)
        else:
            print('Close Figure 1 window to continue...')
        plt.show()
    mpl.rcParams['axes.titlesize'] = current_title_size
    mpl.rcParams['axes.labelsize'] = current_label_size
    mpl.rcParams['xtick.labelsize'] = current_xlabel_size
    mpl.rcParams['ytick.labelsize'] = current_ylabel_size
    return 0


def check_tolerance(expected, actual, abstol=0, reltol=1e-4):
    if abs(expected - actual) > abstol + reltol*abs(expected):
        raise ValueError("|%.3g - %.3g| > %.3g" %
                         (expected, actual, abstol + reltol*abs(expected)))


class pTest(object):
    def __init__(self, filename, record, throw=False):
        self.record = record
        self.data = []
        self.ntest = 0
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
            if self.throw:
                raise ValueError(msg)
            print(msg)
        if self.record:
            self.file.close()

    def check(self, value, abs_tol=0, rel_tol=1e-3):
        if self.record:
            self.file.write('%e\n' % value)
        else:
            if self.ntest >= self.size:
                raise IndexError('no data available for test %d' % self.ntest)
            else:
                expected = self.data[self.ntest]
                try:
                    check_tolerance(expected, value, abs_tol, rel_tol)
                except ValueError as e:
                    self.failed += 1
                    msg = ('+++ test %d failed:' % self.ntest) + str(e)
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


def check_status(handle, stack = None):
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
            repr(line) + ' of ' + file + '\n' + \
            'the reconstruction engine output may provide more information'
        raise error(errorMsg)


def try_calling(returned_handle):
    check_status(returned_handle, inspect.stack()[1])
    pyiutil.deleteDataHandle(returned_handle)


def assert_validity(object, type):
    assert isinstance(object, type)
    assert object.handle is not None


def assert_validities(x, y):
    assert type(x) == type(y)
    assert x.handle is not None
    assert y.handle is not None


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
