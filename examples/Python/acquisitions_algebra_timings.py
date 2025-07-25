import numpy
import sirf
import sirf.STIR as pet
import sirf.Gadgetron as mr
import sys
import timeit

from sirf.Utilities import examples_data_path, existing_filepath
from sirf.SIRF import norm, dot, copyto

narg = len(sys.argv)
if narg < 2:
    print('Usage: acq_algebra_timings <acquisition data filepath> [<number of tests>, <check_norm>]')
    exit()

data_file = sys.argv[1]

if sys.argv[1].endswith('.h5'):
    data_path = examples_data_path('MR')
    data_file = existing_filepath(data_path, sys.argv[1])
    acq = mr.AcquisitionData(data_file)
    mr.AcquisitionData.set_storage_scheme('memory')
    mod = 'mr'
else:
    data_path = examples_data_path('PET')
    data_file = existing_filepath(data_path, sys.argv[1])
    acq = pet.AcquisitionData(data_file)
    pet.AcquisitionData.set_storage_scheme('memory')
    mod = 'pet'

if narg > 2:
    ntests = int(sys.argv[2])
else:
    ntests = 1
check_norm = (ntests == 1)

x = acq + 0
y = x + 0
z = x + 0

if mod == 'mr':
    x_view = mr.AcquisitionDataView(x)
    y_view = mr.AcquisitionDataView(y)
    z_view = mr.AcquisitionDataView(z)
else:
    x_view = x.asarray(copy=False)
    y_view = y.asarray(copy=False)
    z_view = z.asarray(copy=False)

tests = 8
'''
tests:
x * 2
x + 2
x + y
x * y
x / y
norm
sum
dot
'''

view_t = numpy.zeros(tests)
sirf_t = numpy.zeros(tests)


for test in range(ntests):

    start = timeit.default_timer()
    y = 2*x
    elapsed = timeit.default_timer() - start
    sirf_t[0] += elapsed
    norm_y = y.norm()

    start = timeit.default_timer()
    copyto(y_view, x_view)
    y_view *= 2
    elapsed = timeit.default_timer() - start
    view_t[0] += elapsed
    if check_norm:
        print(f'norm(y): {norm_y} {norm(y)} {norm(y_view)}')

    start = timeit.default_timer()
    y = x + 2
    elapsed = timeit.default_timer() - start
    sirf_t[1] += elapsed
    norm_y = y.norm()

    start = timeit.default_timer()
    copyto(y_view, x_view)
    y_view += 2
    elapsed = timeit.default_timer() - start
    view_t[1] += elapsed
    if check_norm:
        print(f'norm(y): {norm_y} {norm(y)} {norm(y_view)}')

    start = timeit.default_timer()
    z = x + y
    elapsed = timeit.default_timer() - start
    sirf_t[2] += elapsed
    norm_z = z.norm()

    start = timeit.default_timer()
    copyto(z_view, x_view)
    z_view += y_view
    elapsed = timeit.default_timer() - start
    view_t[2] += elapsed
    if check_norm:
        print(f'norm(z): {norm_z} {norm(z)} {norm(z_view)}')

    start = timeit.default_timer()
    z = x * y
    elapsed = timeit.default_timer() - start
    sirf_t[3] += elapsed
    norm_z = z.norm()

    start = timeit.default_timer()
    copyto(z_view, x_view)
    z_view *= y_view
    elapsed = timeit.default_timer() - start
    view_t[3] += elapsed
    if check_norm:
        print(f'norm(z): {norm_z} {norm(z)} {norm(z_view)}')

    if mod == 'pet':
        y = y.maximum(1e-20)
        y_view = y.asarray(copy=False)

    start = timeit.default_timer()
    z = x / y
    elapsed = timeit.default_timer() - start
    sirf_t[4] += elapsed
    norm_z = z.norm()

    start = timeit.default_timer()
    copyto(z_view, x_view)
    z_view /= y_view
    elapsed = timeit.default_timer() - start
    view_t[4] += elapsed
    if check_norm:
        print(f'norm(z): {norm_z} {norm(z)} {norm(z_view)}')

    start = timeit.default_timer()
    s = norm(x)
    elapsed = timeit.default_timer() - start
    sirf_t[5] += elapsed

    start = timeit.default_timer()
    t = norm(x_view)
    elapsed = timeit.default_timer() - start
    view_t[5] += elapsed
    if ntests == 1:
        print(f'norm: {s} {t}')

    start = timeit.default_timer()
    s = x.sum()
    elapsed = timeit.default_timer() - start
    sirf_t[6] += elapsed

    start = timeit.default_timer()
    t = x_view.sum()
    elapsed = timeit.default_timer() - start
    view_t[6] += elapsed
    if ntests == 1:
        print(f'sum: {s} {t}')

    start = timeit.default_timer()
    s = x.dot(y)
    elapsed = timeit.default_timer() - start
    sirf_t[7] += elapsed

    start = timeit.default_timer()
    t = dot(x_view, y_view)
    elapsed = timeit.default_timer() - start
    view_t[7] += elapsed
    if ntests == 1:
        print(f'dot: {s} {t}')

sirf_t /= ntests
view_t /= ntests
print('test     sirf    sirf with views')
print(f'x * 2  {sirf_t[0]:.2e}     {view_t[0]:.2e}')
print(f'x + 2  {sirf_t[1]:.2e}     {view_t[1]:.2e}')
print(f'x + y  {sirf_t[2]:.2e}     {view_t[2]:.2e}')
print(f'x * y  {sirf_t[3]:.2e}     {view_t[3]:.2e}')
print(f'x / y  {sirf_t[4]:.2e}     {view_t[4]:.2e}')
print(f'norm   {sirf_t[5]:.2e}     {view_t[5]:.2e}')
print(f'sum    {sirf_t[6]:.2e}     {view_t[6]:.2e}')
print(f'dot    {sirf_t[7]:.2e}     {view_t[7]:.2e}')

