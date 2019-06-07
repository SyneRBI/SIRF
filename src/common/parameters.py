import inspect

import sirf.select_module as select_module
import sirf.pyiutilities as pyiutil
from sirf.Utilities import check_status

exec('from sirf.' + select_module.module + ' import setParameter, parameter')

def _set_parameter(hs, set, par, hv, stack = None):
    if stack is None:
        stack = inspect.stack()[1]
    h = setParameter(hs, set, par, hv)
    check_status(h, stack)
    pyiutil.deleteDataHandle(h)
def set_char_par(handle, set, par, value):
    h = pyiutil.charDataHandle(value)
    _set_parameter(handle, set, par, h, inspect.stack()[1])
    pyiutil.deleteDataHandle(h)
def set_int_par(handle, set, par, value):
    h = pyiutil.intDataHandle(value)
    _set_parameter(handle, set, par, h, inspect.stack()[1])
    pyiutil.deleteDataHandle(h)
def set_float_par(handle, set, par, value):
    h = pyiutil.floatDataHandle(value)
    _set_parameter(handle, set, par, h, inspect.stack()[1])
    pyiutil.deleteDataHandle(h)
def char_par(handle, set, par):
    h = parameter(handle, set, par)
    check_status(h)
    value = pyiutil.charDataFromHandle(h)
    pyiutil.deleteDataHandle(h)
    return value
def int_par(handle, set, par):
    h = parameter(handle, set, par)
    check_status(h, inspect.stack()[1])
    value = pyiutil.intDataFromHandle(h)
    pyiutil.deleteDataHandle(h)
    return value
def int_pars(handle, set, par, n):
    h = parameter(handle, set, par)
    check_status(h)
    for i in range(n):
        value += (pyiutil.intDataItemFromHandle(h, i),)
    pyiutil.deleteDataHandle(h)
    return value
def uint16_pars(handle, set, par, n):
    h = parameter(handle, set, par)
    check_status(h)
    value = ()
    for i in range(n):
        value += (pyiutil.uint16DataItemFromHandle(h, i),)
    pyiutil.deleteDataHandle(h)
    return value
def uint32_pars(handle, set, par, n):
    h = parameter(handle, set, par)
    check_status(h)
    value = ()
    for i in range(n):
        value += (pyiutil.uint32DataItemFromHandle(h, i),)
    pyiutil.deleteDataHandle(h)
    return value
def uint64_pars(handle, set, par, n):
    h = parameter(handle, set, par)
    check_status(h)
    value = ()
    for i in range(n):
        value += (pyiutil.uint64DataItemFromHandle(h, i),)
    pyiutil.deleteDataHandle(h)
    return value
def float_par(handle, set, par):
    h = parameter(handle, set, par)
    check_status(h)
    v = pyiutil.floatDataFromHandle(h)
    pyiutil.deleteDataHandle(h)
    return v
def float_pars(handle, set, par, n):
    h = parameter(handle, set, par)
    check_status(h)
    value = ()
    for i in range(n):
        value += (pyiutil.floatDataItemFromHandle(h, i),)
    pyiutil.deleteDataHandle(h)
    return value
def parameter_handle(hs, set, par):
    handle = parameter(hs, set, par)
    check_status(handle, inspect.stack()[1])
    return handle
