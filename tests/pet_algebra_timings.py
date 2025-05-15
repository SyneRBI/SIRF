import numpy
import sirf.STIR as pet
import sys
import timeit

from sirf.Utilities import examples_data_path, existing_filepath

narg = len(sys.argv)
if narg < 2:
    print('Usage: pet_algebra_timings <data filepath> [<number of tests>]')
    exit()

data_path = examples_data_path('PET')
data_file = existing_filepath(data_path, sys.argv[1])

argv1 = sys.argv[1]
if argv1.endswith('.hv'):
    x = pet.ImageData(sys.argv[1])
else:
    pet.AcquisitionData.set_storage_scheme('memory')
    x = pet.AcquisitionData(sys.argv[1])

# provide direct access to data
y = x.clone()
z = x.clone()
xa = x.asarray()
ya = y.asarray()
za = z.asarray()

if narg > 2:
    ntests = int(sys.argv[2])
else:
    ntests = 1

tests = 5
'''
tests:
x * 2
x + 2
x + y
x * y
x / y
'''

sarr_t = numpy.zeros(tests)
sirf_t = numpy.zeros(tests)
asar_t = numpy.zeros(tests)
nump_t = numpy.zeros(tests)
fill_t = numpy.zeros(tests)

for test in range(ntests):

    start = timeit.default_timer()
    numpy.copyto(ya, xa)
    ya *= 2
    elapsed = timeit.default_timer() - start
    sarr_t[0] += elapsed
    start = timeit.default_timer()
    y = 2*x
    elapsed = timeit.default_timer() - start
    sirf_t[0] += elapsed
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

    start = timeit.default_timer()
    numpy.copyto(ya, xa)
    ya += 2
    elapsed = timeit.default_timer() - start
    sarr_t[1] += elapsed
    start = timeit.default_timer()
    y = x + 2
    elapsed = timeit.default_timer() - start
    sirf_t[1] += elapsed
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

    start = timeit.default_timer()
    numpy.copyto(za, xa)
    za += ya
    elapsed = timeit.default_timer() - start
    sarr_t[2] += elapsed
    start = timeit.default_timer()
    xy = x + y
    elapsed = timeit.default_timer() - start
    sirf_t[2] += elapsed
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
    xy.fill(z_arr)
    elapsed = timeit.default_timer() - start
    fill_t[2] += elapsed

    start = timeit.default_timer()
    numpy.copyto(za, xa)
    za *= ya
    elapsed = timeit.default_timer() - start
    sarr_t[3] += elapsed
    start = timeit.default_timer()
    xy = x * y
    elapsed = timeit.default_timer() - start
    sirf_t[3] += elapsed
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
    xy.fill(z_arr)
    elapsed = timeit.default_timer() - start
    fill_t[3] += elapsed

    y += 1e-20
    start = timeit.default_timer()
    numpy.copyto(za, xa)
    za /= ya
    elapsed = timeit.default_timer() - start
    sarr_t[4] += elapsed
    y = y.maximum(1e-20)
    start = timeit.default_timer()
    xy = x / y
    elapsed = timeit.default_timer() - start
    sirf_t[4] += elapsed
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
    xy.fill(z_arr)
    elapsed = timeit.default_timer() - start
    fill_t[4] += elapsed

print('test     sirf    as_array +   numpy  +   fill   =  total   with asarray')
numpy_tot = asar_t[0] + nump_t[0] + fill_t[0]
print(f'x * 2  {sirf_t[0]:.2e}  {asar_t[0]:.2e} + {nump_t[0]:.2e} + {fill_t[0]:.2e} = {numpy_tot:.2e}   {sarr_t[0]:.2e}')
numpy_tot = asar_t[1] + nump_t[1] + fill_t[1]
print(f'x + 2  {sirf_t[1]:.2e}  {asar_t[1]:.2e} + {nump_t[1]:.2e} + {fill_t[1]:.2e} = {numpy_tot:.2e}   {sarr_t[1]:.2e}')
numpy_tot = asar_t[2] + nump_t[2] + fill_t[2]
print(f'x + y  {sirf_t[2]:.2e}  {asar_t[2]:.2e} + {nump_t[2]:.2e} + {fill_t[2]:.2e} = {numpy_tot:.2e}   {sarr_t[2]:.2e}')
numpy_tot = asar_t[3] + nump_t[3] + fill_t[3]
print(f'x * y  {sirf_t[3]:.2e}  {asar_t[3]:.2e} + {nump_t[3]:.2e} + {fill_t[3]:.2e} = {numpy_tot:.2e}   {sarr_t[3]:.2e}')
numpy_tot = asar_t[4] + nump_t[4] + fill_t[4]
print(f'x / y  {sirf_t[4]:.2e}  {asar_t[4]:.2e} + {nump_t[4]:.2e} + {fill_t[4]:.2e} = {numpy_tot:.2e}   {sarr_t[4]:.2e}')

