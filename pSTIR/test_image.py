import numpy
import pylab
import stir

try:
    printer = stir.printerTo('stdout')

    voxel_dim = (100, 100, 50)
    voxel_size = (2, 2, 1)

    image = stir.Image()
    image.initialise(voxel_dim, voxel_size)

    shape = stir.EllipsoidalCylinder()
    shape.set_length(40)
    shape.set_radii((30, 10))
    shape.set_origin((0, 20, 10))
    image.add_shape(shape, scale = 1)

    shape.set_radii((10, 10))
    shape.set_origin((20, -10, 10))
    image.add_shape(shape, scale = 2)

    shape.set_origin((-20, -10, 10))
    image.add_shape(shape, scale = 0.25)

    data = image.density()

    nz = data.shape[0]
    ny = data.shape[1]
    nx = data.shape[2]
    print(nx, ny, nz)

    while True:
        s = str(input('enter z-coordinate: '))
        if len(s) < 1:
            break
        z = int(s)
        if z < 0 or z >= nz:
            break
        pylab.figure(z)
        pylab.imshow(data[z,:,:])
        pylab.show()

except stir.error as err:
    print('STIR exception occured:\n', err.value)
