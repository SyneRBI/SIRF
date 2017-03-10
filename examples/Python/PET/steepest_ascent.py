'''Steepest ascent demo

Usage:
  steepest_ascent [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: my_forward_projection.hs]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/PET
                              subfolder of SIRF root folder
  -t <step>, --tau=<step>     steepest ascent step size parameter,
                              use a negative value to opt for the optimal value
                              [default: -1]
  -s <nstp>, --steps=<nstp>   number of steepest descent steps [default: 3]
  -v, --verbose               verbose
  -e <engn>, --engine=<engn>  reconstruction engine [default: Stir]
'''

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

import scipy.optimize

# import engine module
exec('from p' + args['--engine'] + ' import *')

# process command-line options
tau = float(args['--tau'])
steps = int(args['--steps'])
verbose = args['--verbose']
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

    # direct all information printing to a file
    info_printer = printerTo('info.txt', INFO_CHANNEL)
    # direct all warning printing to a file
    warning_printer = printerTo('warn.txt', WARNING_CHANNEL)
    # direct all error printing to stdout
    error_printer = printerTo('stdout', ERROR_CHANNEL)

    # create matrix to be used by the acquisition model
    matrix = RayTracingMatrix().set_num_tangential_LORs(2)

    # create acquisition model
    am = AcquisitionModelUsingMatrix(matrix)

    # PET acquisition data to be read from the file specified by --file option
    print('raw data: %s' % raw_data_file)
    ad = AcquisitionData(raw_data_file)

    # create filter
    filter = CylindricFilter()

    # create initial image estimate
    nx = 111
    ny = 111
    nz = 31
    image_size = (nx, ny, nz)
    voxel_size = (3, 3, 3.375)
    image = ImageData()
    image.initialise(image_size, voxel_size)
    image.fill(1.0)
    filter.apply(image)

    # create objective function
    obj_fun = PoissonLogLh_LinModMean_AcqMod()
    obj_fun.set_acquisition_model(am)
    obj_fun.set_acquisition_data(ad)
    obj_fun.set_num_subsets(12)
    obj_fun.set_up(image)

    # plot the initial image
    idata = image.as_array()
    show(1, 'initial image', idata[20,:,:])

    print('computing initial objective function value...')
    print('objective function value: %e' % (obj_fun.value(image)))

    if verbose:
        disp = 3
        print('NOTE: below f(x) is the negative of the objective function value')
    else:
        disp = 0
    eps = 1e-6

    for iter in range(steps):

        # obtain gradient
        grad = obj_fun.gradient(image, 0)
        filter.apply(grad)
        gdata = grad.as_array()

        # find maximal steepest descent step parameter t in image + t*grad 
        # such that the new image remains positive
        max_image = idata.max()
        max_grad = abs(gdata).max()
        gdata[abs(gdata) < eps] = eps
        r = idata/gdata
        gdata[idata <= 0] = 0
        d = r[r < 0]
        if d.shape[0] > 0:
            maxstep = -d.max()
        else:
            maxstep = abs(tau)*max_image/max_grad

        if tau < 0:
            # find the optimal t
            fun = lambda t: -obj_fun.value(image.fill(idata + t*gdata))
            t = scipy.optimize.fminbound \
                (fun, 0, maxstep, xtol = 1e-4, maxfun = 3, disp = disp)
        else:
            # t is such that the relative change in image is not greater than tau
            t = tau*max_image/max_grad
            if t > maxstep:
                t = maxstep

        # perform steepest descent step
        print('step %d, max change in image %e' % (iter, t*max_grad))
        idata = idata + t*gdata

        # filter the new image
        image.fill(idata)
        filter.apply(image)
        idata = image.as_array()

        # plot the new image
        show(iter + 1, 'current image', idata[20,:,:])

        # quit if the new image has negative values
        min_image = idata.min()
        if min_image < -eps:
            print('image minimum is negative: %e' % min_image)
            break

    if tau > 0 or disp == 0:
        print('computing final objective function value...')
        print('objective function value: %e' % (obj_fun.value(image)))

# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
    print('done')
except error as err:
    # display error information
    print('STIR exception occured: %s' % err.value)
