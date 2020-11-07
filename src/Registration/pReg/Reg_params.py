import inspect

import sirf.select_module as select_module
import sirf.pyiutilities as pyiutil
from sirf.Utilities import check_status

from sirf.pyreg import setParameter, parameter


def set_parameter(hs, group, par, hv, stack = None):
    if stack is None:
        stack = inspect.stack()[1]
    h = setParameter(hs, group, par, hv)
    check_status(h, stack)
    pyiutil.deleteDataHandle(h)


def set_bool_par(handle, group, par, value):
    h = pyiutil.charDataHandle(value)
    set_parameter(handle, group, par, h, inspect.stack()[1])
    pyiutil.deleteDataHandle(h)


def set_char_par(handle, group, par, value):
    h = pyiutil.charDataHandle(value)
    set_parameter(handle, group, par, h, inspect.stack()[1])
    pyiutil.deleteDataHandle(h)


def set_int_par(handle, group, par, value):
    h = pyiutil.intDataHandle(value)
    set_parameter(handle, group, par, h, inspect.stack()[1])
    pyiutil.deleteDataHandle(h)


def set_float_par(handle, group, par, value):
    h = pyiutil.floatDataHandle(value)
    set_parameter(handle, group, par, h, inspect.stack()[1])
    pyiutil.deleteDataHandle(h)


def bool_par(handle, group, par):
    h = parameter(handle, group, par)
    check_status(h, inspect.stack()[1])
    value = pyiutil.boolDataFromHandle(h)
    pyiutil.deleteDataHandle(h)
    return value


def char_par(handle, group, par):
    h = parameter(handle, group, par)
    check_status(h)
    value = pyiutil.charDataFromHandle(h)
    pyiutil.deleteDataHandle(h)
    return value


def int_par(handle, group, par):
    h = parameter(handle, group, par)
    check_status(h, inspect.stack()[1])
    value = pyiutil.intDataFromHandle(h)
    pyiutil.deleteDataHandle(h)
    return value


def int_pars(handle, group, par, n):
    h = parameter(handle, group, par)
    check_status(h)
    value = ()
    for i in range(n):
        value += (pyiutil.intDataItemFromHandle(h, i),)
    pyiutil.deleteDataHandle(h)
    return value


def uint16_pars(handle, group, par, n):
    h = parameter(handle, group, par)
    check_status(h)
    value = ()
    for i in range(n):
        value += (pyiutil.uint16DataItemFromHandle(h, i),)
    pyiutil.deleteDataHandle(h)
    return value


def uint32_pars(handle, group, par, n):
    h = parameter(handle, group, par)
    check_status(h)
    value = ()
    for i in range(n):
        value += (pyiutil.uint32DataItemFromHandle(h, i),)
    pyiutil.deleteDataHandle(h)
    return value


def uint64_pars(handle, group, par, n):
    h = parameter(handle, group, par)
    check_status(h)
    value = ()
    for i in range(n):
        value += (pyiutil.uint64DataItemFromHandle(h, i),)
    pyiutil.deleteDataHandle(h)
    return value


def float_par(handle, group, par):
    h = parameter(handle, group, par)
    check_status(h)
    v = pyiutil.floatDataFromHandle(h)
    pyiutil.deleteDataHandle(h)
    return v


def float_pars(handle, group, par, n):
    h = parameter(handle, group, par)
    check_status(h)
    value = ()
    for i in range(n):
        value += (pyiutil.floatDataItemFromHandle(h, i),)
    pyiutil.deleteDataHandle(h)
    return value


def parameter_handle(hs, group, par):
    handle = parameter(hs, group, par)
    check_status(handle, inspect.stack()[1])
    return handle
