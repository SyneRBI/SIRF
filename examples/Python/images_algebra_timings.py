import numpy
import sirf.STIR as pet
import sirf.Gadgetron as mr
import sirf.Reg as reg
import sys
import timeit

from sirf.Utilities import examples_data_path, existing_filepath
from sirf.SIRF import norm, dot, copyto

narg = len(sys.argv)
if narg < 2:
    print('Usage: images_algebra_timings <data filepath> [<number of tests>]')
    exit()

if sys.argv[1].endswith('.h5'):
    data_path = examples_data_path('MR')
    data_file = existing_filepath(data_path, sys.argv[1])
    mod = 'mr'
    x = mr.ImageData(data_file)
    data_type = x.data_type(0)
    if data_type < 5:
        print('integer data not supported, converting to real...')
        x = x.real()
elif sys.argv[1].endswith('.nii') or sys.argv[1].endswith('.nii.gz'):
    data_path = examples_data_path('Registration')
    data_file = existing_filepath(data_path, sys.argv[1])
    mod = 'reg'
    x = reg.ImageData(data_file)
else:
    data_path = examples_data_path('PET')
    data_file = existing_filepath(data_path, sys.argv[1])
    mod = 'pet'
    x = pet.ImageData(data_file)

if narg > 2:
    ntests = int(sys.argv[2])
else:
    ntests = 1
check_norm = narg > 3

if mod == 'pet':
    y = x + 0
    z = x + 0
else:
    y = x.clone()
    z = x.clone()
if mod == 'mr':
    x_view = mr.ImageDataView(x)
    y_view = mr.ImageDataView(y)
    z_view = mr.ImageDataView(z)
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
    norm_y = norm(y)
#    norm_y = y.norm()

    start = timeit.default_timer()
    '''
    we cannot do y_view = x_view * 2 because then y_view would no longer be a view of y,
    and y would remain unchanged - so we employ this workaround:
    '''
    copyto(y_view, x_view)
    y_view *= 2
    elapsed = timeit.default_timer() - start
    view_t[0] += elapsed
    if ntests == 1:
        print(f'norm(y): {norm_y} {norm(y)} {norm(y_view)}')

    start = timeit.default_timer()
    y = x + 2
    elapsed = timeit.default_timer() - start
    sirf_t[1] += elapsed
    norm_y = norm(y)

    start = timeit.default_timer()
    copyto(y_view, x_view)
    y_view += 2
    elapsed = timeit.default_timer() - start
    view_t[1] += elapsed
    if ntests == 1:
        print(f'norm(y): {norm_y} {norm(y)} {norm(y_view)}')

    start = timeit.default_timer()
    z = x + y
    elapsed = timeit.default_timer() - start
    sirf_t[2] += elapsed
    norm_z = norm(z)

    start = timeit.default_timer()
    copyto(z_view, x_view)
    z_view += y_view
    elapsed = timeit.default_timer() - start
    view_t[2] += elapsed
    if ntests == 1:
        print(f'norm(z): {norm_z} {norm(z)} {norm(z_view)}')

    start = timeit.default_timer()
    z = x * y
    elapsed = timeit.default_timer() - start
    sirf_t[3] += elapsed
    norm_z = norm(z)

    start = timeit.default_timer()
    copyto(z_view, x_view)
    z_view *= y_view
    elapsed = timeit.default_timer() - start
    view_t[3] += elapsed
    if ntests == 1:
        print(f'norm(z): {norm_z} {norm(z)} {norm(z_view)}')

    if mod == 'pet':
        y = y.maximum(1e-20)
    start = timeit.default_timer()
    z = x / y
    elapsed = timeit.default_timer() - start
    sirf_t[4] += elapsed
    norm_z = norm(z)

    start = timeit.default_timer()
    copyto(z_view, x_view)
    z_view /= y_view
    elapsed = timeit.default_timer() - start
    view_t[4] += elapsed
    if ntests == 1:
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

