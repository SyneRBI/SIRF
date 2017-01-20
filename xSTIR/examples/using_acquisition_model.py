'''Forward projection demo: creates an image, forward-projects it to simulate
acquisition data and displays

Usage:
  using_acquisition_model [--help | options]

Options:
  -e <engn>, --engine=<engn>  reconstruction engine [default: Stir]
  -p <path>, --path=<path>    sub-path to engine module [default: /xSTIR/pSTIR]

There is an interactive demo with much more documentation on this process.
You probably want to check that instead.
'''

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

import os
import sys
sys.path.append(os.environ.get('SRC_PATH') + args['--path'])
exec('from p' + args['--engine'] + ' import *')

def main():

    # output goes to files
    printer = Printer('info.txt', 'warn.txt', 'errr.txt')

    # create an empty image
    image = Image()
    image_size = (111, 111, 31)
    voxel_size = (3, 3, 3.375)
    image.initialise(image_size, voxel_size)

    # create a shape
    shape = EllipsoidalCylinder()

    # add a shape
    shape.set_length(400)
    shape.set_radii((100, 40))
    shape.set_origin((0, 60, 10))
    image.add_shape(shape, scale = 1)

    # add another shape
    shape.set_radii((30, 30))
    shape.set_origin((60, -30, 10))
    image.add_shape(shape, scale = 1.5)

    # add another shape
    shape.set_origin((-60, -30, 10))
    image.add_shape(shape, scale = 0.75)

    # z-pixel coordinate of the xy-crossection to plot
    z = int(image_size[2]/2)

    # plot the phantom image to be reconstructed
    data = image.as_array()
    pylab.figure(1000)
    pylab.imshow(data[z,:,:])
    print('Figure 1000: exact image - close window to continue')
    pylab.show()

    # define the matrix to be used by the acquisition model
    matrix = RayTracingMatrix()
    matrix.set_num_tangential_LORs(2)

    # define the acquisition model
    am = AcquisitionModelUsingMatrix()
    am.set_matrix(matrix)

    # define a filter
    filter = CylindricFilter()

    # create an initial image estimate
    reconstructedImage = Image()
    reconstructedImage.initialise(image_size, voxel_size)
    reconstructedImage.fill(1.0)
    # apply filter to get a cylindric initial image
    filter.apply(reconstructedImage)

    # plot the initial image
    data = reconstructedImage.as_array()
    pylab.figure(1)
    pylab.imshow(data[z,:,:])
    print('Figure 1: initial image - close window to continue')
    pylab.show()

    print('projecting image...')
    # forward-project the image to obtain 'raw data'
    # 'Utahscat600k_ca_seg4.hs' is used as a template
    templ = AcquisitionData('Utahscat600k_ca_seg4.hs')
    am.set_up(templ, image)
    ad = am.forward(image)
    # if the raw data is very large, it can be stored in a file
    # ad = am.forward(image, 'demo4data.hs')

    print('back-projecting image...')
    # backward-project the computed forward projection
    update = am.backward(ad)


try:
    main()
except error as err:
    print('exception occured: %s' % err.value)
