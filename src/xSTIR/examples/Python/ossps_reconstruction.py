'''OSSPS reconstruction demo with user-controlled iterations

Usage:
  ossps_reconstruction [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: Utahscat600k_ca_seg4.hs]
  -P <path>, --path=<path>    path to data files, defaults to data/examples/PET
                              subfolder of SIRF root folder
  -g <file>, --init=<file>    initial image guess file
                              [default: test_image_PM_QP_6.hv]
  -p <fact>, --penf=<fact>    penalty factor [default: 0]
  -s <subs>, --subs=<subs>    number of subsets [default: 4]
  -i <iter>, --iter=<iter>    number of iterations [default: 2]
  -e <engn>, --engine=<engn>  reconstruction engine [default: Stir]
'''

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

pen_factor = args['--penf']
num_subsets = int(args['--subs'])
num_subiterations = int(args['--iter'])

# import engine module
exec('from p' + args['--engine'] + ' import *')

def main():

    # direct all printing to this file
    printer = printerTo('output.txt')

    # select acquisition model that implements the geometric
    # forward projection by a ray tracing matrix multiplication
    am = AcquisitionModelUsingMatrix\
         (RayTracingMatrix().set_num_tangential_LORs(2))

    # locate the input data file folder
    data_path = args['--path']
    if data_path is None:
        data_path = pet_data_path()

    # PET acquisition data to be read from the file specified by --file option
    raw_data_file = existing_filepath(data_path, args['--file'])
    ad = AcquisitionData(raw_data_file)

    # define objective function to be maximized as
    # Poisson logarithmic likelihood with linear model for mean
    obj_fun = PoissonLogLh_LinModMean_AcqMod()
    obj_fun.set_acquisition_model(am)
    obj_fun.set_acquisition_data(ad)
    fact = float(pen_factor)
    obj_fun.set_prior(QuadraticPrior().set_penalisation_factor(fact))

    # select Ordered Subsets Separable Paraboloidal Surrogate
    # as the reconstruction algorithm
    recon = OSSPSReconstruction()
    recon.set_num_subsets(num_subsets)
    recon.set_num_subiterations(num_subiterations)
    recon.set_objective_function(obj_fun)

    # read an initial estimate for the reconstructed image from the file
    # specified by --init option
    init_image_file = existing_filepath(data_path, args['--init'])
    image = ImageData(init_image_file)

    # set up the reconstructor
    print('setting up, please wait...')
    recon.set_up(image)

    # in order to see the reconstructed image evolution
    # take over the control of the iterative process
    # rather than allow recon.reconstruct to do all job at once
    for iter in range(num_subiterations):
        print('\n------------- Subiteration %d' % recon.get_subiteration_num())
        # perform an iteration
        recon.update(image)
        # plot the current image
        data = image.as_array()
        pylab.figure(iter + 1)
        pylab.imshow(data[10,:,:])
        print('close Figure %d window to continue' % iter)
        pylab.show()

    # let the user inspect any z-crossections of the image they want to
    data = image.as_array()
    nz = data.shape[0]
    print('Enter z-coordinate of the slice to view it')
    print('(a value outside the range [1 : %d] will stop the loop)'%nz)
    while True:
        s = str(input('z-coordinate: '))
        if len(s) < 1:
            break
        z = int(s)
        if z < 1 or z > nz:
            break
        pylab.figure(z)
        pylab.imshow(data[z - 1,:,:])
        print('close Figure %d window to continue' % z)
        pylab.show()

# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
except error as err:
    # display error information
    print('STIR exception occured: %s' % err.value)
