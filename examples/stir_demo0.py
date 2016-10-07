import argparse
import numpy
import pylab
import os
import sys
sys.path.append(os.environ.get('CSTIR_SRC') + '/../pSTIR')

from pStir import *

parser = argparse.ArgumentParser(description = \
'''
Phantom creating demo
''')
args = parser.parse_args()

def main():

    image = PETImage()
    image_dim = (100, 100, 50)
    voxel_size = (0.8, 0.8, 1)
    image.initialise(image_dim, voxel_size)
    data = image.as_array()
    pylab.figure(1)
    pylab.imshow(data[0,:,:])
    print('close Figure 1 window to continue')
    pylab.show()

    shape = EllipsoidalCylinder()

    shape.set_length(40)
    shape.set_radii((30, 10))
    shape.set_origin((0, 20, 10))
    print('adding shape 1...')
    image.add_shape(shape, scale = 1)
    data = image.as_array()
    pylab.figure(2)
    pylab.imshow(data[0,:,:])
    print('close Figure 2 window to continue')
    pylab.show()

    shape.set_radii((10, 10))
    shape.set_origin((20, -10, 10))
    print('adding shape 2...')
    image.add_shape(shape, scale = 2)
    data = image.as_array()
    pylab.figure(3)
    pylab.imshow(data[0,:,:])
    print('close Figure 3 window to continue')
    pylab.show()

    shape.set_origin((-20, -10, 10))
    print('adding shape 3...')
    image.add_shape(shape, scale = 0.25)
    data = image.as_array()
    pylab.figure(4)
    pylab.imshow(data[0,:,:])
    print('close Figure 4 window to continue')
    pylab.show()

    print('removing shape 3...')
    image.add_shape(shape, scale = -0.25)
    data = image.as_array()
    pylab.figure(5)
    pylab.imshow(data[0,:,:])
    print('close Figure 5 window to continue')
    pylab.show()

try:
    main()

except error as err:
    print('STIR exception occured: %s' % err.value)
