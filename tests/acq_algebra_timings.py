import numpy
import sirf.STIR as pet
import sirf.Gadgetron as mr
import sys
import timeit

from sirf.Utilities import examples_data_path, existing_filepath

narg = len(sys.argv)
if narg < 2:
    print('Usage: acq_algebra_timings <acquisition data filepath> [<number of tests>, <check_norm>]')
    exit()

pet.AcquisitionData.set_storage_scheme('memory')

if sys.argv[1].endswith('.h5'):
    data_path = examples_data_path('MR')
    data_file = existing_filepath(data_path, sys.argv[1])
    acq = mr.AcquisitionData(data_file)
else:
    data_path = examples_data_path('PET')
    data_file = existing_filepath(data_path, sys.argv[1])
    acq = pet.AcquisitionData(data_file)

if narg > 2:
    ntests = int(sys.argv[2])
else:
    ntests = 1
check_norm = narg > 3

x = acq.clone()

#print(f'x.norm() = {x.norm()}')

tests = 5
'''
tests:
x * 2
x + 2
x + y
x * y
x / y
'''

sirf_t = numpy.zeros(tests)
asar_t = numpy.zeros(tests)
nump_t = numpy.zeros(tests)
fill_t = numpy.zeros(tests)

for test in range(ntests):

    start = timeit.default_timer()
    y = 2*x
    elapsed = timeit.default_timer() - start
    sirf_t[0] += elapsed
    norm_y = y.norm()
    start = timeit.default_timer()
    x_arr = x.as_array()
    elapsed = timeit.default_timer() - start
    asar_t[0] += elapsed
    start = timeit.default_timer()
    y_arr = 2*x_arr
    elapsed = timeit.default_timer() - start
    nump_t[0] += elapsed
    start = timeit.default_timer()
    y.fill(y_arr)
    elapsed = timeit.default_timer() - start
    fill_t[0] += elapsed
    if check_norm:
        print(f'norm(y): {norm_y} {y.norm()}')

    start = timeit.default_timer()
    y = x + 2
    elapsed = timeit.default_timer() - start
    sirf_t[1] += elapsed
    norm_y = y.norm()
    start = timeit.default_timer()
    x_arr = x.as_array()
    elapsed = timeit.default_timer() - start
    asar_t[1] += elapsed
    start = timeit.default_timer()
    y_arr = x_arr + 2
    elapsed = timeit.default_timer() - start
    nump_t[1] += elapsed
    start = timeit.default_timer()
    y.fill(y_arr)
    elapsed = timeit.default_timer() - start
    fill_t[1] += elapsed
    if check_norm:
        print(f'norm(y): {norm_y} {y.norm()}')

    start = timeit.default_timer()
    z = x + y
    elapsed = timeit.default_timer() - start
    sirf_t[2] += elapsed
    norm_z = z.norm()
    start = timeit.default_timer()
    x_arr = x.as_array()
    y_arr = y.as_array()
    elapsed = timeit.default_timer() - start
    asar_t[2] += elapsed
    start = timeit.default_timer()
    z_arr = x_arr + y_arr
    elapsed = timeit.default_timer() - start
    nump_t[2] += elapsed
    start = timeit.default_timer()
    z.fill(z_arr)
    elapsed = timeit.default_timer() - start
    fill_t[2] += elapsed
    if check_norm:
        print(f'norm(z): {norm_z} {z.norm()}')

    start = timeit.default_timer()
    z = x * y
    elapsed = timeit.default_timer() - start
    sirf_t[3] += elapsed
    norm_z = z.norm()
    start = timeit.default_timer()
    x_arr = x.as_array()
    y_arr = y.as_array()
    elapsed = timeit.default_timer() - start
    asar_t[3] += elapsed
    start = timeit.default_timer()
    z_arr = x_arr * y_arr
    elapsed = timeit.default_timer() - start
    nump_t[3] += elapsed
    start = timeit.default_timer()
    z.fill(z_arr)
    elapsed = timeit.default_timer() - start
    fill_t[3] += elapsed
    if check_norm:
        print(f'norm(z): {norm_z} {z.norm()}')

    start = timeit.default_timer()
    z = x / y
    elapsed = timeit.default_timer() - start
    sirf_t[4] += elapsed
    norm_z = z.norm()
    start = timeit.default_timer()
    x_arr = x.as_array()
    y_arr = y.as_array()
    elapsed = timeit.default_timer() - start
    asar_t[4] += elapsed
    start = timeit.default_timer()
    z_arr = x_arr / y_arr
    elapsed = timeit.default_timer() - start
    nump_t[4] += elapsed
    start = timeit.default_timer()
    z.fill(z_arr)
    elapsed = timeit.default_timer() - start
    fill_t[4] += elapsed
    if check_norm:
        print(f'norm(z): {norm_z} {z.norm()}')

print('test     sirf    as_array +   numpy  +   fill   =  total')
numpy_tot = asar_t[0] + nump_t[0] + fill_t[0]
print(f'x * 2  {sirf_t[0]:.2e}  {asar_t[0]:.2e} + {nump_t[0]:.2e} + {fill_t[0]:.2e} = {numpy_tot:.2e}')
numpy_tot = asar_t[1] + nump_t[1] + fill_t[1]
print(f'x + 2  {sirf_t[1]:.2e}  {asar_t[1]:.2e} + {nump_t[1]:.2e} + {fill_t[1]:.2e} = {numpy_tot:.2e}')
numpy_tot = asar_t[2] + nump_t[2] + fill_t[2]
print(f'x + y  {sirf_t[2]:.2e}  {asar_t[2]:.2e} + {nump_t[2]:.2e} + {fill_t[2]:.2e} = {numpy_tot:.2e}')
numpy_tot = asar_t[3] + nump_t[3] + fill_t[3]
print(f'x * y  {sirf_t[3]:.2e}  {asar_t[3]:.2e} + {nump_t[3]:.2e} + {fill_t[3]:.2e} = {numpy_tot:.2e}')
numpy_tot = asar_t[4] + nump_t[4] + fill_t[4]
print(f'x / y  {sirf_t[4]:.2e}  {asar_t[4]:.2e} + {nump_t[4]:.2e} + {fill_t[4]:.2e} = {numpy_tot:.2e}')

