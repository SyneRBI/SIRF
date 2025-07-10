import numpy
import sirf.STIR as pet
import sirf.Gadgetron as mr
import sirf.Reg as reg
import sys
import timeit

from sirf.Utilities import examples_data_path, existing_filepath

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
elif sys.argv[1].endswith('.nii'):
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

tests = 5
'''
tests:
x * 2
x + 2
x + y
x * y
x / y
'''

view_t = numpy.zeros(tests)
sirf_t = numpy.zeros(tests)


def copy_view(mod, x, y):
    if mod == 'mr':
        y.copy(x)
    else:
        numpy.copyto(y, x)


for test in range(ntests):

    start = timeit.default_timer()
    y = 2*x
    elapsed = timeit.default_timer() - start
    sirf_t[0] += elapsed
    norm_y = y.norm()

    start = timeit.default_timer()
    copy_view(mod, x_view, y_view)
    y_view *= 2
    elapsed = timeit.default_timer() - start
    view_t[0] += elapsed
    if ntests == 1:
        print(f'norm(y): {norm_y} {y.norm()}')

    start = timeit.default_timer()
    y = x + 2
    elapsed = timeit.default_timer() - start
    sirf_t[1] += elapsed
    norm_y = y.norm()

    start = timeit.default_timer()
    copy_view(mod, x_view, y_view)
    y_view += 2
    elapsed = timeit.default_timer() - start
    view_t[1] += elapsed
    if ntests == 1:
        print(f'norm(y): {norm_y} {y.norm()}')

    start = timeit.default_timer()
    z = x + y
    elapsed = timeit.default_timer() - start
    sirf_t[2] += elapsed
    norm_z = z.norm()

    start = timeit.default_timer()
    copy_view(mod, x_view, z_view)
    z_view += y_view
    elapsed = timeit.default_timer() - start
    view_t[2] += elapsed
    if ntests == 1:
        print(f'norm(z): {norm_z} {z.norm()}')

    start = timeit.default_timer()
    z = x * y
    elapsed = timeit.default_timer() - start
    sirf_t[3] += elapsed
    norm_z = z.norm()

    start = timeit.default_timer()
    copy_view(mod, x_view, z_view)
    z_view *= y_view
    elapsed = timeit.default_timer() - start
    view_t[3] += elapsed
    if ntests == 1:
        print(f'norm(z): {norm_z} {z.norm()}')

    if mod == 'pet':
        y = y.maximum(1e-20)
    start = timeit.default_timer()
    z = x / y
    elapsed = timeit.default_timer() - start
    sirf_t[4] += elapsed
    norm_z = z.norm()

    start = timeit.default_timer()
    copy_view(mod, x_view, z_view)
    elapsed = timeit.default_timer() - start
    view_t[4] += elapsed
    if ntests == 1:
        print(f'norm(z): {norm_z} {z.norm()}')

print('test     sirf    sirf with views')
print(f'x * 2  {sirf_t[0]:.2e}     {view_t[0]:.2e}')
print(f'x + 2  {sirf_t[1]:.2e}     {view_t[1]:.2e}')
print(f'x + y  {sirf_t[2]:.2e}     {view_t[2]:.2e}')
print(f'x * y  {sirf_t[3]:.2e}     {view_t[3]:.2e}')
print(f'x / y  {sirf_t[4]:.2e}     {view_t[4]:.2e}')

