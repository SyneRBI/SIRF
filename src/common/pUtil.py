'''Utilities used by all engines
'''

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
## Copyright 2015 - 2017 University College London.
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

import matplotlib.pyplot as plt
import os
import pyiutil

class error(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return '??? ' + repr(self.value)

def petmr_data_path(petmr):
    data_path = '/data/examples/' + petmr.upper()
    SIRF_PATH = os.environ.get('SIRF_PATH')
    if SIRF_PATH is not None:
        return SIRF_PATH + data_path
    SRC_PATH = os.environ.get('SRC_PATH')
    if SRC_PATH is None:
        errorMsg = 'Path to raw data files not found'
        raise error(errorMsg)
    return SRC_PATH + '/SIRF' + data_path

def existing_filepath(data_path, file_name):
    full_name = data_path + '/' + file_name
    if not os.path.isfile(full_name):
        raise error('file %s not found' % full_name)
    return full_name

def check_status(handle):
    if pyiutil.executionStatus(handle) != 0:
        msg = pyiutil.executionError(handle)
        file = pyiutil.executionErrorFile(handle)
        line = pyiutil.executionErrorLine(handle)
        errorMsg = \
            repr(msg) + ' exception thrown at line ' + \
            repr(line) + ' of ' + file
        raise error(errorMsg)

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

def show_3D_array\
    (array, tile_shape = None, scale = None, \
     suptitle = None, titles = None, label = None, show = True):
    import math
    import numpy
    if tile_shape is None:
        nz = array.shape[0]
        ny = array.shape[1]
        nx = array.shape[2]
        rows = int(round(math.sqrt(nz*nx/ny)))
        if rows < 1:
            rows = 1
        if rows > nz:
            rows = nz
        cols = (nz - 1)//rows + 1
    else:
        rows, cols = tile_shape
        assert rows*cols >= array.shape[0],\
                "tile rows x columns must equal the 3rd dim extent of array"
    if scale is None:
        vmin = numpy.amin(array)
        vmax = numpy.amax(array)
    else:
        vmin, vmax = scale
    fig = plt.figure()
    if suptitle is not None:
        fig.suptitle(suptitle, fontsize = 16)
    for z in range(array.shape[0]):
        ax = fig.add_subplot(rows, cols, z+1)
        if titles is None:
            if label is None:
                ax.set_title('%d' % (z + 1))
            else:
                ax.set_title(label + (' %d' % (z + 1)), fontsize = 8)
        else:
            ax.set_title(titles[z])
        ax.set_axis_off()
        imgplot = ax.imshow(array[z,:,:], vmin=vmin, vmax=vmax)
    print('close figure 1 to continue')
    if show:
        plt.show()

