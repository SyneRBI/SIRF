# -*- coding: utf-8 -*-

"""
sirf.Gadgetron Views Test.
Compares the performance of two approaches that avoid copying data between C++ and Python:
- C++ algebra that is optimised by using templated loops
- Python algebra that directly accesses C++ data via NumPy Array API views

Usage:
  mr_timings [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: simulated_MR_2D_cartesian.h5]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/MR
                              subfolder of SIRF root folder
  -t <tsts>, --tests=<tsts>    number of tests [default: 1]
  --non-interactive           do not show plots
"""

import numpy
import timeit
import importlib

from sirf.SIRF import dot, copyto
from sirf.Utilities import examples_data_path, existing_filepath

__version__ = "0.1.0"
from docopt import docopt
args = docopt(__doc__, version=__version__)
#print(args)

# process command-line options
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = examples_data_path('MR')
input_file = existing_filepath(data_path, data_file)
ntests = int(args['--tests'])

# import engine module
mr = importlib.import_module('sirf.Gadgetron')

# obtain acquisitions and image data
acq = mr.AcquisitionData(input_file)
processed_data = mr.preprocess_acquisition_data(acq)
recon = mr.FullySampledReconstructor()
recon.set_input(processed_data)
print('---\n reconstructing...')
recon.process()
img = recon.get_output()

tests = 5
'''
tests:
x * 2
x + 2
x + y
x * y
x / y
'''

print('\n=== Comparing acquisitions algebra norms and timings...\n')

ax = acq
ay = ax + 0
az = ay + 0
vx = mr.AcquisitionDataView(ax)
vy = mr.AcquisitionDataView(ay)
vz = mr.AcquisitionDataView(az)

temp_t = numpy.zeros(tests)
view_t = numpy.zeros(tests)

for test in range(ntests):

    start = timeit.default_timer()
    ay = 2*ax
    elapsed = timeit.default_timer() - start
    temp_t[0] += elapsed
    norm_y = ay.norm()
    start = timeit.default_timer()
    vy.copy(vx)
    vy *= 2
    elapsed = timeit.default_timer() - start
    view_t[0] += elapsed
    if test == 0:
        print(f'norm(x * 2): {norm_y} {vy.norm()}')

    start = timeit.default_timer()
    ay = ax + 2
    elapsed = timeit.default_timer() - start
    temp_t[1] += elapsed
    norm_y = ay.norm()
    start = timeit.default_timer()
    vy.copy(vx)
    vy += 2
    elapsed = timeit.default_timer() - start
    view_t[1] += elapsed
    if test == 0:
        print(f'norm(x + 2): {norm_y} {vy.norm()}')

    start = timeit.default_timer()
    az = ax + ay
    elapsed = timeit.default_timer() - start
    temp_t[2] += elapsed
    norm_z = az.norm()
    start = timeit.default_timer()
    vz.copy(vx)
    vz += vy
    elapsed = timeit.default_timer() - start
    view_t[2] += elapsed
    if test == 0:
        print(f'norm(x + y): {norm_z} {vz.norm()}')

    start = timeit.default_timer()
    az = ax * ay
    elapsed = timeit.default_timer() - start
    temp_t[3] += elapsed
    norm_z = az.norm()
    start = timeit.default_timer()
    vz.copy(vx)
    vz *= vy
    elapsed = timeit.default_timer() - start
    view_t[3] += elapsed
    if test == 0:
        print(f'norm(x * y): {norm_z} {vz.norm()}')

    ay += 1e-10
    start = timeit.default_timer()
    az = ax / ay
    elapsed = timeit.default_timer() - start
    temp_t[4] += elapsed
    norm_z = az.norm()
    start = timeit.default_timer()
    vz.copy(vx)
    vz /= vy
    elapsed = timeit.default_timer() - start
    view_t[4] += elapsed
    if test == 0:
        print(f'norm(x / y): {norm_z} {vz.norm()}')

view_t /= ntests
temp_t /= ntests

print('\ntest    templates    views')
print(f'x * 2   {temp_t[0]:.2e}    {view_t[0]:.2e}')
print(f'x + 2   {temp_t[1]:.2e}    {view_t[1]:.2e}')
print(f'x + y   {temp_t[2]:.2e}    {view_t[2]:.2e}')
print(f'x * y   {temp_t[3]:.2e}    {view_t[3]:.2e}')
print(f'x / y   {temp_t[4]:.2e}    {view_t[4]:.2e}')

print('\n=== Comparing images algebra norms and timings...\n')

ix = img
iy = ix + 0
iz = iy + 0
vx = mr.ImageDataView(ix)
vy = mr.ImageDataView(iy)
vz = mr.ImageDataView(iz)

temp_t = numpy.zeros(tests)
view_t = numpy.zeros(tests)

for test in range(ntests):

    start = timeit.default_timer()
    iy = 2*ix
    elapsed = timeit.default_timer() - start
    temp_t[0] += elapsed
    norm_y = iy.norm()
    start = timeit.default_timer()
    vy.copy(vx)
    vy *= 2
    elapsed = timeit.default_timer() - start
    view_t[0] += elapsed
    if test == 0:
        print(f'norm(x * 2): {norm_y} {vy.norm()}')

    start = timeit.default_timer()
    iy = ix + 2
    elapsed = timeit.default_timer() - start
    temp_t[1] += elapsed
    norm_y = iy.norm()
    start = timeit.default_timer()
    vy.copy(vx)
    vy += 2
    elapsed = timeit.default_timer() - start
    view_t[1] += elapsed
    if test == 0:
        print(f'norm(x + 2): {norm_y} {vy.norm()}')

    start = timeit.default_timer()
    iz = ix + iy
    elapsed = timeit.default_timer() - start
    temp_t[2] += elapsed
    norm_z = iz.norm()
    start = timeit.default_timer()
    vz.copy(vx)
    vz += vy
    elapsed = timeit.default_timer() - start
    view_t[2] += elapsed
    if test == 0:
        print(f'norm(x + y): {norm_z} {vz.norm()}')

    start = timeit.default_timer()
    iz = ix * iy
    elapsed = timeit.default_timer() - start
    temp_t[3] += elapsed
    norm_z = iz.norm()
    start = timeit.default_timer()
    vz.copy(vx)
    vz *= vy
    elapsed = timeit.default_timer() - start
    view_t[3] += elapsed
    if test == 0:
        print(f'norm(x * y): {norm_z} {vz.norm()}')

    ay += 1e-10
    start = timeit.default_timer()
    iz = ix / iy
    elapsed = timeit.default_timer() - start
    temp_t[4] += elapsed
    norm_z = iz.norm()
    start = timeit.default_timer()
    vz.copy(vx)
    vz /= vy
    elapsed = timeit.default_timer() - start
    view_t[4] += elapsed
    if test == 0:
        print(f'norm(x / y): {norm_z} {vz.norm()}')

view_t /= ntests
temp_t /= ntests

print('\ntest    templates    views')
print(f'x * 2   {temp_t[0]:.2e}    {view_t[0]:.2e}')
print(f'x + 2   {temp_t[1]:.2e}    {view_t[1]:.2e}')
print(f'x + y   {temp_t[2]:.2e}    {view_t[2]:.2e}')
print(f'x * y   {temp_t[3]:.2e}    {view_t[3]:.2e}')
print(f'x / y   {temp_t[4]:.2e}    {view_t[4]:.2e}')

print('\n=== done with %s' % __file__)
