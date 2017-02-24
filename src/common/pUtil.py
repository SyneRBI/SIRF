'''Utilities used by all engines
'''
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

