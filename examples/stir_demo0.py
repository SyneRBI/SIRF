import numpy
try:
    import pylab
    HAVE_PYLAB = True
except:
    HAVE_PYLAB = False
import sys
sys.path.append('..\pSTIR')
import stir

try:
    image = stir.Image()
    image_dim = (100, 100, 50)
    voxel_size = (0.8, 0.8, 1)
    image.initialise(image_dim, voxel_size)
    if HAVE_PYLAB:
        data = image.as_array()
        pylab.imshow(data[0,:,:])
        pylab.show()

    shape = stir.EllipsoidalCylinder()

    shape.set_length(40)
    shape.set_radii((30, 10))
    shape.set_origin((0, 20, 10))
    print('adding shape 1...')
    image.add_shape(shape, scale = 1)
    if HAVE_PYLAB:
        data = image.as_array()
        pylab.imshow(data[0,:,:])
        pylab.show()

    shape.set_radii((10, 10))
    shape.set_origin((20, -10, 10))
    print('adding shape 2...')
    image.add_shape(shape, scale = 2)
    if HAVE_PYLAB:
        data = image.as_array()
        pylab.imshow(data[0,:,:])
        pylab.show()

    shape.set_origin((-20, -10, 10))
    print('adding shape 3...')
    image.add_shape(shape, scale = 0.25)
    if HAVE_PYLAB:
        data = image.as_array()
        pylab.imshow(data[0,:,:])
        pylab.show()

    print('removing shape 3...')
    image.add_shape(shape, scale = -0.25)
    if HAVE_PYLAB:
        data = image.as_array()
        pylab.imshow(data[0,:,:])
        pylab.show()

except stir.error as err:
    print('STIR exception occured:\n', err.value)
