# import graphics library
import pylab
# pylab draws numpy arrays
import numpy
# import STIR interface library
import sys
sys.path.append('..\pSTIR')
import stir

# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    # direct all diagnostic printing to a file
    printer = stir.printerTo('stir_demo3.txt')

    # create matrix to be used by projectors
    matrix = stir.RayTracingMatrix()
    matrix.set_num_tangential_LORs(2)
    n = matrix.get_num_tangential_LORs()

    # create acquisition model
    am = stir.AcquisitionModelUsingMatrix()
    am.set_matrix(matrix)

    # read acquisition model data
    ad = stir.AcquisitionData('Utahscat600k_ca_seg4.hs')
##    ad.read_from_file('Utahscat600k_ca_seg4.hs')

    # create prior
    prior = stir.QuadraticPrior()
    prior.set_penalisation_factor(0.5)

    # create objective function
    obj_fun = stir.PoissonLogLh_LinModMean_AcqModData()
    obj_fun.set_sensitivity_filename('RPTsens_seg3_PM.hv')
    obj_fun.set_recompute_sensitivity(True)
    obj_fun.set_use_subset_sensitivities(False)
    obj_fun.set_zero_seg0_end_planes(True)
    obj_fun.set_max_segment_num_to_process(3)
    obj_fun.set_acquisition_model(am)
    obj_fun.set_acquisition_data(ad)
    obj_fun.set_prior(prior)

    # create OSSPS reconstructor
    recon = stir.OSSPSReconstruction()
    recon.set_output_filename_prefix('reconstructedImage')
    recon.set_num_subsets(4)
    recon.set_num_subiterations(8)
    recon.set_save_interval(8)
    recon.set_relaxation_parameter(2)
    recon.set_objective_function(obj_fun)

    # read an initial estimate for the reconstructed image from a file
    image = stir.Image()
    image.read_from_file('test_image_PM_QP_6.hv')

    # set up the reconstructor
    recon.set_up(image)

    # obtain the subiterations range
    start = recon.get_start_subiteration_num()
    stop = recon.get_num_subiterations()
    print('subiteration range:', start, stop)

    # in order to see the reconstructed image evolution
    # take over the control of the iterative process
    # rather than allow recon.reconstruct to do all job at once

    # not needed (set by set_up), just to show we can start at any iteration 
    recon.set_subiteration_num(start)
    for iter in range(start, stop + 1):
        print('\n--------------------- Subiteration ',\
              recon.get_subiteration_num())
        # perform an iteration
        recon.update(image)
        # plot the current image
        data = image.as_array()
        pylab.figure(iter)
        pylab.imshow(data[10,:,:])
        pylab.show()

    # compare the reconstructed image to the expected image
    expectedImage = stir.Image()
    expectedImage.read_from_file('test_image_OSSPS_PM_QP_8.hv')
    diff = expectedImage.diff_from(image)
    print('difference from expected image:', diff)

    # let the user inspect any z-crossections of the image they want to
    data = image.as_array()
    nz = data.shape[0]
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
    # display error information
    print('STIR exception occured:\n', err.value)
