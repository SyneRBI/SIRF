import argparse
import numpy
import os
try:
    import pylab
    HAVE_PYLAB = True
except:
    HAVE_PYLAB = False
import sys
sys.path.append(os.environ.get('CSTIR_SRC') + '/../pSTIR')
import stir
import time

parser = argparse.ArgumentParser(description = \
'''
OSMAPOSL reconstruction demo with all parameters defined in the script
and user-controlled iterations
''')
parser.add_argument('-n', '--normalisation', action = 'store', type = float, \
    default = 1, help = 'normalisation, default 1')
parser.add_argument('-a', '--additive', action = 'store', type = float, \
    default = 0, help = 'additive term, default 1')
args = parser.parse_args()

def main():

    # direct all information printing to a file
    #info_printer = stir.printerTo('stdout', stir.INFO_CHANNEL)
    info_printer = stir.printerTo('stir_demo2info.txt', stir.INFO_CHANNEL)
    # direct all warning printing to a file
    warning_printer = stir.printerTo('stir_demo2warn.txt', stir.WARNING_CHANNEL)
    # direct all error printing to stdout
    error_printer = stir.printerTo('stdout', stir.ERROR_CHANNEL)

    exact_image = stir.Image('my_image.hv')
    image = exact_image.get_empty_copy()
    image.fill(1.0)

    acq_templ = stir.AcquisitionData('my_forward_projection.hs')

    # create matrix to be used by the acquisition model
    matrix = stir.RayTracingMatrix()
    matrix.set_num_tangential_LORs(2)

    # create acquisition model
    #am = stir.AcquisitionModelUsingMatrix()
    am = stir.PETAcquisitionModelUsingMatrix()
    am.set_matrix(matrix)
    am.set_up(acq_templ, exact_image)

    print('projecting image...')
    ad = am.forward(exact_image)

    # define a filter
    filter = stir.CylindricFilter()

    # define the objective function
    obj_fun = stir.PoissonLogLh_LinModMean_AcqModData()
    #obj_fun.set_max_segment_num_to_process(3)
    obj_fun.set_pet_acquisition_model(am)
    obj_fun.set_acquisition_data(ad)

    num_subiterations = 2

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
    print('setting up reconstructor, please wait...')
    recon.set_up(image)

    for iter in range(1, num_subiterations + 1):
        print('\n--------------------- Subiteration %d'\
              % recon.get_subiteration_num())
        # perform an iteration
        recon.update(image)

    expected_image = image.clone()
    image.fill(1.0)

    add = stir.AcquisitionData(acq_templ)
    bkg = stir.AcquisitionData(acq_templ)
    nrm = stir.AcquisitionData(acq_templ)
    add.fill(args.additive)
    bkg.fill(args.additive)
    nrm.fill(args.normalisation)

    print('\n--- testing normalisation only...')
    am.set_normalisation(nrm)
    am.set_up(acq_templ, exact_image)
    print('projecting image...')
    new_ad = am.forward(exact_image)
    obj_fun.set_pet_acquisition_model(am)
    obj_fun.set_acquisition_data(new_ad)

    print('setting up reconstructor, please wait...')
    recon.set_up(image)

    for iter in range(1, num_subiterations + 1):
        print('\n--------------------- Subiteration %d'\
              % recon.get_subiteration_num())
        # perform an iteration
        recon.update(image)

    # compare the reconstructed image to the expected image
    diff = expected_image.diff_from(image)
    print('\n--- difference from expected image: %e' % diff)

    print('\n--- testing normalisation and additive term...')
    am.set_additive_term(add)
    am.set_normalisation(nrm)
    am.set_up(acq_templ, exact_image)
    print('projecting image...')
    new_ad = am.forward(exact_image)
    obj_fun.set_pet_acquisition_model(am)
    obj_fun.set_acquisition_data(new_ad)

    image = exact_image.clone()
    print('setting up reconstructor, please wait...')
    recon.set_up(image)

    for iter in range(1, num_subiterations + 1):
        print('\n--------------------- Subiteration %d'\
              % recon.get_subiteration_num())
        # perform an iteration
        recon.update(image)

    # compare the reconstructed image to the expected image
    diff = exact_image.diff_from(image)
    print('\n--- difference from expected image: %e' % diff)

## not possible to test - background term not implemented in STIR yet
##    print('\n--- testing background term...')
##    add.fill(0.0)
##    nrm.fill(1.0)
##    am.set_background_term(bkg)
##    am.set_normalisation(nrm)
##    am.set_up(acq_templ, exact_image)
##    print('projecting image...')
##    new_ad = am.forward(exact_image)
##    obj_fun.set_pet_acquisition_model(am)
##    obj_fun.set_acquisition_data(new_ad)
##
##    image = exact_image.clone()
##    print('setting up reconstructor, please wait...')
##    recon.set_up(image)
##
##    for iter in range(1, num_subiterations + 1):
##        print('\n--------------------- Subiteration %d'\
##              % recon.get_subiteration_num())
##        # perform an iteration
##        recon.update(image)
##
##    # compare the reconstructed image to the expected image
##    diff = exact_image.diff_from(image)
##    print('\n--- difference from expected image: %e' % diff)

# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
except stir.error as err:
    # display error information
    print('STIR exception occured: %s\n' % err.value)
