'''Forward projection demo: creates an image, forward-projects it to simulate
acquisition data and displays

Usage:
  using_acquisition_model [--help | options]

Options:
  -f <file>, --file=<file>    raw data file [default: Utahscat600k_ca_seg4.hs]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/PET
                              subfolder of SIRF root folder
  -e <engn>, --engine=<engn>  reconstruction engine [default: Stir]

There is an interactive demo with much more documentation on this process.
You probably want to check that instead.
'''

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

# import engine module
exec('from p' + args['--engine'] + ' import *')

# process command-line options
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = petmr_data_path('pet')
raw_data_file = existing_filepath(data_path, data_file)

def main():

    # output goes to files
    printer = Printer('info.txt', 'warn.txt', 'errr.txt')

    # create an empty image
    image = ImageData()
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
    reconstructedImage = ImageData()
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
    # raw_data_file is used as a template
    templ = AcquisitionData(raw_data_file)
    am.set_up(templ, image)
    ad = am.forward(image)
    # if the raw data is very large, it can be stored in a file
    # ad = am.forward(image, 'proj_data.hs')

    print('back-projecting the forward projection...')
    # backward-project the computed forward projection
    update = am.backward(ad)

try:
    main()
except error as err:
    print('exception occured: %s' % err.value)
