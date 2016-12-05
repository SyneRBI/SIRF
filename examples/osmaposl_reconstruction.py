"""OSMAPOSL reconstruction demo with user-controlled iterations

Usage:
  osmaposl_reconstruction [--help | options]

Options:
  -e=<e>, --engine=<e>  reconstruction engine [default: Stir]
  --log=<lvl>  CRITICAL|ERROR|WARN(ING)|[default: INFO]|DEBUG
"""
# @comment 2016-12-02 Casper dC-L <imaging@caspersci.uk.to>
# argparse replaced with docopt
import logging
# @comment 2016-12-02 Casper dC-L <imaging@caspersci.uk.to>
# add built-in logging framework
import matplotlib.pyplot as plt
# @comment 2016-12-02 Casper dC-L <imaging@caspersci.uk.to>
# replace pylab with plt <- matplotlib.pyplot
__version__ = "0.1.0"
# @comment 2016-12-02 Casper dC-L <imaging@caspersci.uk.to>
# version number added


# a simplistic example of user's involvement in the reconstruction
def my_image_data_processor(image_array, im_num):
    log = logging.getLogger(__name__)
    # plot the current estimate of the image
    plt.figure(im_num)
    plt.title('image estimate %d' % im_num)
    plt.imshow(image_array[20, :, :])
    log.info('close Figure %d window to continue' % im_num)
    plt.show()
    # image is not modified in this simplistic example - but might have been
    return image_array


def run(args):
    exec("from p" + args["--engine"] + " import *")
    # @comment 2016-12-02 Casper dC-L <imaging@caspersci.uk.to>
    # The following is preferable but we may not want it for simplicity
    #   from importlib import import_module
    #   pEngine = import_module("p" + args["--engine"])
    # Then use:
    #   pEngine.AcquisitionModelUsingMatrix, etc...

    log = logging.getLogger(__name__)

    # select acquisition model that implements the geometric
    # forward projection by a ray tracing matrix multiplication
    am = AcquisitionModelUsingMatrix(RayTracingMatrix())

    # indicate that acquisition data source is this file
    # (which contains forward projection of the actual image
    # stored in file my_image.hv - see below)
    ad = AcquisitionData('my_forward_projection.hs')

    # create initial image estimate compatible with acquisition data
    # and initialize each voxel to 1.0
    image = ad.create_empty_image(1.0)

    # define objective function to be maximized as
    # Poisson logarithmic likelihood with linear model for mean
    # that uses acquisition model for computing the gradient
    obj_fun = PoissonLogLh_LinModMean_AcqMod()
    obj_fun.set_acquisition_model(am)
    obj_fun.set_acquisition_data(ad)
    # @comment 2016-12-02 Casper dC-L <imaging@caspersci.uk.to>
    # do not like apparently thoroughly inconsistent combination of:
    # - title case
    # - lowercase
    # - underscores

    # select Ordered Subsets Maximum A-Posteriori One Step Late
    # as the reconstruction algorithm
    recon = OSMAPOSLReconstruction()
    recon.set_objective_function(obj_fun)
    recon.set_num_subsets(12)

    # set up the reconstructor
    log.info('setting up, please wait...')
    recon.set_up(image)

    # in order to see the reconstructed image evolution
    # take over the control of the iterative process
    # rather than allow recon.reconstruct to do all job at once
    recon.set_current_estimate(image)
    for iteration in range(2):
        log.debug('\n------------- iteration %d' % iteration)
        # perform an iteration
        recon.update_current_estimate()
        # copy current image estimate into python array to inspect/process
        image_array = recon.get_current_estimate().as_array()
        # apply your image data processor/visualizer
        my_image_array = my_image_data_processor(image_array, iteration + 1)
        # fill the current image estimate with new data
        recon.get_current_estimate().fill(my_image_array)

    # compare the reconstructed image to the actual image
    actual_image_array = Image('my_image.hv').as_array()
    plt.figure(iteration + 1)
    plt.title('reconstructed image')
    plt.imshow(image_array[20, :, :])
    plt.figure(0)
    plt.title('actual image')
    plt.imshow(actual_image_array[20, :, :])
    log.info('close Figure windows %d and %d to continue' % (
        0, iteration + 1))
    plt.show()


def main():
    from docopt import docopt
    args = docopt(__doc__, version=__version__)
    logging.basicConfig(level=getattr(logging, args["--log"], logging.INFO))
    # @comment 2016-12-02 Casper dC-L <imaging@caspersci.uk.to>
    # print replaced with logging. customisations replaced with:
    # Optional: direct all information printing and warnings to files
    # (only need to do this once in sirf base package).
    #   log = logging.getLogger(__name__)
    #   infoHand = logging.FileHandler('info.txt')
    #   infoHand.setLevel(logging.INFO)
    #   log.addHandler(infoHand)
    #   warnHand = logging.FileHandler('warn.txt')
    #   warnHand.setLevel(logging.WARN)
    #   log.addHandler(warnHand)

    run(args)


# if anything goes wrong, an exception will be thrown
# (cf. Error Handling section in the spec)
if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        # display error information
        from sys import stderr
        # @comment 2016-12-02 Casper dC-L <imaging@caspersci.uk.to>
        # error to console error stream
        stderr.write(__file__ + ':ERROR:' + str(e))
