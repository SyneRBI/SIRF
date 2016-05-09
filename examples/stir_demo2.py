import numpy
try:
    import pylab
    HAVE_PYLAB = True
except:
    HAVE_PYLAB = False
import sys
sys.path.append('../pSTIR')
import stir
import time

# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    # direct all information printing to a file
    info_printer = stir.printerTo('stir_demo2info.txt', stir.INFO_CHANNEL)
    # direct all warning printing to a file
    warning_printer = stir.printerTo('stir_demo2warn.txt', stir.WARNING_CHANNEL)
    # direct all error printing to stdout
    error_printer = stir.printerTo('stdout', stir.ERROR_CHANNEL)

    # create matrix to be used by the acquisition model
    matrix = stir.RayTracingMatrix()
    matrix.set_num_tangential_LORs(2)

    # create acquisition model
    am = stir.AcquisitionModelUsingMatrix()
    am.set_matrix(matrix)

    # read acquisition model data
    ad = stir.AcquisitionData('my_forward_projection.hs')
    #amd.read_from_file('my_forward_projection.hs')

    # create prior
    prior = stir.QuadraticPrior()
    prior.set_penalisation_factor(0.001)

    # create filter
    filter = stir.CylindricFilter()

    # create initial image estimate
    image_size = (111, 111, 31)
    voxel_size = (3, 3, 3.375)
    image = stir.Image()
    image.initialise(image_size, voxel_size)
    image.fill(1.0)
##    filter.set_strictly_less_than_radius(False)
##    filter.apply(image)
##    filter.set_strictly_less_than_radius(True)

    # create objective function
    obj_fun = stir.PoissonLogLh_LinModMean_AcqModData()
    obj_fun.set_zero_seg0_end_planes(True)
    obj_fun.set_max_segment_num_to_process(3)
    obj_fun.set_acquisition_model(am)
    obj_fun.set_acquisition_data(ad)
    obj_fun.set_prior(prior)

    num_subiterations = 6

    # create OSMAPOSL reconstructor
    recon = stir.OSMAPOSLReconstruction()
    recon.set_objective_function(obj_fun)
    recon.set_MAP_model('multiplicative')
    recon.set_num_subsets(12)
    recon.set_num_subiterations(num_subiterations)
    recon.set_save_interval(num_subiterations)
    recon.set_inter_iteration_filter_interval(1)
    recon.set_inter_iteration_filter(filter)
    recon.set_output_filename_prefix('reconstructedImage')

    # set up the reconstructor
    print('setting up, please wait...')
    recon.set_up(image)

    if HAVE_PYLAB:
        # plot the initial image
        data = image.as_array()
        pylab.figure(1)
        pylab.imshow(data[20,:,:])
        print('Figure 1: initial image - close window to continue')
        pylab.show()

    #recon.reconstruct(image)

    # in order to see the reconstructed image evolution
    # take over the control of the iterative process
    # rather than allow recon.reconstruct to do all job at once
    for iter in range(1, num_subiterations + 1):
        print('\n------------- Subiteration %d' % recon.get_subiteration_num())
        # perform an iteration
        recon.update(image)
        if HAVE_PYLAB:
            # plot the current image
            data = image.as_array()
            pylab.figure(iter + 1)
            pylab.imshow(data[20,:,:])
            print('close Figure %d window to continue' % (iter + 1))
            pylab.show()
        # image can be post-processed
        #filter.apply(image)

    # compare the reconstructed image to the expected image
    expectedImage = stir.Image('expected_image.hv')
    diff = expectedImage.diff_from(image)
    print('difference from expected image: %e' % diff)

    # compare the reconstructed image to the exact image
    exactImage = stir.Image('my_image.hv')
    x_data = exactImage.as_array()
    data = image.as_array()

    if HAVE_PYLAB:
        pylab.figure(1000)
        pylab.imshow(data[20,:,:])
        pylab.figure(1001)
        pylab.imshow(x_data[20,:,:])
        print('close Figure 1001 window, then')
        print('close Figure 1000 window to continue')
        pylab.show()

except stir.error as err:
    # display error information
    print('STIR exception occured:\n', err.value)
