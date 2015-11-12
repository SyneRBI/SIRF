import numpy
import pylab
import stir

try:
    printer = stir.printerTo('stdout')

    image0 = stir.Image()
    voxel_dim = (100, 100, 50)
    voxel_size = (2, 2, 1)
    image0.initialise(voxel_dim, voxel_size)

    shape = stir.EllipsoidalCylinder()

    shape.set_length(400)
    shape.set_radii((400, 400))
    shape.set_origin((150, 0, 0))

    shape.set_length(40)
    shape.set_radii((30, 10))
    shape.set_origin((0, 20, 10))
    image0.add_shape(shape, scale = 1)

    shape.set_radii((10, 10))
    shape.set_origin((20, -10, 10))
    image0.add_shape(shape, scale = 2)

    shape.set_origin((-20, -10, 10))
    image0.add_shape(shape, scale = 0.25)

    data0 = image0.density()
    nz = data0.shape[0]
    while True:
        s = str(input('enter z-coordinate: '))
        if len(s) < 1:
            break
        z = int(s)
        if z < 0 or z >= nz:
            break
        pylab.figure(z)
        pylab.imshow(data0[z,:,:])
        pylab.show()

    image = stir.Image()
    voxel_dim = (60, 60, 31)
    voxel_size = (4.44114, 4.44114, 3.375)
    voxels = stir.Voxels(voxel_dim, voxel_size)
##    image.initialise(voxel_dim, voxel_size)
    image.initialise(60, 60, 31, 4.44114, 4.44114, 3.375)
    image.fill(1.0)

##    filter = stir.TruncateToCylindricalFOVImageProcessor()
    filter = stir.CylindricFilter()
    filter.set_strictly_less_than_radius(False)
    filter.apply(image)

    image_x = stir.Image()
    image_x.read_from_file('my_uniform_image_circular.hv')

    image_c = image.clone()
    diff = image_x.diff_from(image_c)
    print('difference from cloned image:', diff)

    image_e = image.get_empty_copy(2)

    data = image_e.density()
    print(data.min(), data.max())

    data_x = image_x.density()

    nz = data.shape[0]
    ny = data.shape[1]
    nx = data.shape[2]
    print(nx, ny, nz)

    diff = image_x.diff_from(image_c)
    print('difference from expected image:', diff)

    while True:
        s = str(input('enter z-coordinate: '))
        if len(s) < 1:
            break
        z = int(s)
        if z < 0 or z >= nz:
            break
        pylab.figure(z)
        pylab.imshow(data[z,:,:])
        pylab.figure(1000000 + z)
        pylab.imshow(data_x[z,:,:])
        pylab.show()

except stir.error as err:
    print('STIR exception occured:\n', err.value)
