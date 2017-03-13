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

def show(fig, title, data):
    pylab.figure(fig)
    pylab.title(title)
    pylab.imshow(data)
    pylab.colorbar()
    print('close window to continue')
    pylab.show()

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

    # add a shape to the image
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

    # plot the phantom image
    data = image.as_array()
    show(1, 'Figure 1: phantom image', data[z,:,:])

    # define the acquisition model
    am = AcquisitionModelUsingMatrix()

    print('projecting image...')
    # forward-project the image to obtain simulated acquisition data
    # data from raw_data_file is used as a template
    templ = AcquisitionData(raw_data_file)
    am.set_up(templ, image)
    ad = am.forward(image)
    # if the projection data is very large, it can be stored in a file
    # ad = am.forward(image, 'proj_data.hs')

    # plot simulated acquisition data
    adata = ad.as_array()
    show(2, 'Figure 2: forward projection', adata[z,:,:])

    print('back-projecting the forward projection...')
    # backward-project the computed forward projection
    update = am.backward(ad)

    data = update.as_array()
    show(3, 'Figure 3: back projection', data[z,:,:])

try:
    main()
    print('done')
except error as err:
    print('exception occured: %s' % err.value)
