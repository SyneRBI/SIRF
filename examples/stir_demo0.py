import argparse
import numpy
import os
try:
    import pylab
    HAVE_PYLAB = True
except:
    HAVE_PYLAB = False
import sys
sys.path.append(os.environ.get('CSTIR') + '/../pSTIR')
import stir

parser = argparse.ArgumentParser(description = \
'''
Phantom creating demo
''')
args = parser.parse_args()

def main():

    image = stir.Image()
    image_dim = (100, 100, 50)
    voxel_size = (0.8, 0.8, 1)
    image.initialise(image_dim, voxel_size)
    if HAVE_PYLAB:
        data = image.as_array()
        pylab.figure(1)
        pylab.imshow(data[0,:,:])
        print('close Figure 1 window to continue')
        pylab.show()

    shape = stir.EllipsoidalCylinder()

    shape.set_length(40)
    shape.set_radii((30, 10))
    shape.set_origin((0, 20, 10))
    print('adding shape 1...')
    image.add_shape(shape, scale = 1)
    if HAVE_PYLAB:
        data = image.as_array()
        pylab.figure(2)
        pylab.imshow(data[0,:,:])
        print('close Figure 2 window to continue')
        pylab.show()

    shape.set_radii((10, 10))
    shape.set_origin((20, -10, 10))
    print('adding shape 2...')
    image.add_shape(shape, scale = 2)
    if HAVE_PYLAB:
        data = image.as_array()
        pylab.figure(3)
        pylab.imshow(data[0,:,:])
        print('close Figure 3 window to continue')
        pylab.show()

    shape.set_origin((-20, -10, 10))
    print('adding shape 3...')
    image.add_shape(shape, scale = 0.25)
    if HAVE_PYLAB:
        data = image.as_array()
        pylab.figure(4)
        pylab.imshow(data[0,:,:])
        print('close Figure 4 window to continue')
        pylab.show()

    print('removing shape 3...')
    image.add_shape(shape, scale = -0.25)
    if HAVE_PYLAB:
        data = image.as_array()
        pylab.figure(5)
        pylab.imshow(data[0,:,:])
        print('close Figure 5 window to continue')
        pylab.show()

try:
    main()

except stir.error as err:
    print('STIR exception occured:\n', err.value)
