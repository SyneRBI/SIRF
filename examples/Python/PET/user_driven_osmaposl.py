'''User-driven OSMAPOSL reconstruction

Usage:
  user_driven_osmaposl [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: my_forward_projection.hs]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/PET
                              subfolder of SIRF root folder
  -s <subs>, --subs=<subs>    number of subsets [default: 12]
  -i <iter>, --iter=<iter>    number of iterations [default: 2]
  -e <engn>, --engine=<engn>  reconstruction engine [default: Stir]
'''

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

# import engine module
exec('from p' + args['--engine'] + ' import *')

# process command-line options
num_subsets = int(args['--subs'])
num_subiterations = int(args['--iter'])
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = petmr_data_path('pet')
raw_data_file = existing_filepath(data_path, data_file)

def main():

    # direct all information and error printing to files
    info_printer = printerTo('info.txt', INFO_CHANNEL)
    warning_printer = printerTo('warn.txt', WARNING_CHANNEL)
    # direct all error printing to stdout
    error_printer = printerTo('stdout', ERROR_CHANNEL)

    # create acquisition model
    am = AcquisitionModelUsingMatrix()

    # PET acquisition data to be read from the file specified by --file option
    print('raw data: %s' % raw_data_file)
    ad = AcquisitionData(raw_data_file)

    # create filter
    filter = CylindricFilter()

    # create initial image estimate
    image_size = (111, 111, 31)
    voxel_size = (3, 3, 3.375)
    image = ImageData()
    image.initialise(image_size, voxel_size)
    image.fill(1.0)

    # create prior
    prior = QuadraticPrior()
    prior.set_penalisation_factor(0.5)

    # create objective function
    obj_fun = PoissonLogLh_LinModMean_AcqMod()
    obj_fun.set_acquisition_model(am)
    obj_fun.set_acquisition_data(ad)
    obj_fun.set_num_subsets(num_subsets)
    obj_fun.set_up(image)

    for iter in range(1, num_subiterations + 1):
        print('\n------------- Subiteration %d' % iter) 

        # select subset
        subset = iter - 1

        # get sensitivity as ImageData
        ss_img = obj_fun.get_subset_sensitivity(subset)

        # get gradient (without penalty) + sensitivity as ImageData
        # (back projection of the ratio of measured to estimated acquisition data)
        grad_img = obj_fun.get_gradient_plus_sensitivity_no_penalty(image, subset)

        # get gradient of prior as ImageData
        pgrad_img = prior.get_gradient(image)

        # copy to Python arrays
        image_arr = image.as_array()
        ss_arr = ss_img.as_array()
        grad_arr = grad_img.as_array()
        pgrad_arr = pgrad_img.as_array()

        # update image data
        ss_arr[ss_arr < 1e-6] = 1e-6 # avoid division by zero
        update = grad_arr/(ss_arr + pgrad_arr/num_subsets)
        image_arr = image_arr*update

        # fill current image with new values
        image.fill(image_arr)

        # apply filter
        filter.apply(image)

        # show current image at z = 20
        image_arr = image.as_array()
        pylab.figure(iter)
        pylab.title('Image at z = 20, iteration %d' % iter)
        pylab.imshow(image_arr[20,:,:])
        print('close Figure %d window to continue' % iter)
        pylab.show()

# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
except error as err:
    # display error information
    print('STIR exception occured: %s' % err.value)
