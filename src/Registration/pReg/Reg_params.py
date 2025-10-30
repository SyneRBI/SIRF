"""Internal module for passing/getting Python parameters to/from C/C++."""
from sirf.Utilities import Param
from sirf.pyreg import setParameter, parameter

par = Param(setParameter, parameter)
# backward compatibility
set_parameter = par.set
set_bool_par = par.set_bool
set_char_par = par.set_char
set_int_par = par.set_int
set_float_par = par.set_float
#set_double_par = par.set_double
bool_par = par.get_bool
char_par = par.get_char
int_par = par.get_int
size_t_par = par.get_size_t
int_pars = par.get_ints
uint16_pars = par.get_uint16s
uint32_pars = par.get_uint32s
uint64_pars = par.get_uint64s
float_par = par.get_float
float_pars = par.get_floats
#double_par = par.get_double
parameter_handle = par.get_handle
