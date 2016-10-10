import argparse
import numpy
import pylab
import os
import sys
sys.path.append(os.environ.get('CSTIR_SRC') + '/../pSTIR')
import time

from pStir import *

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
    info_printer = printerTo('info.txt', INFO_CHANNEL)
    # direct all warning printing to a file
    warning_printer = printerTo('warn.txt', WARNING_CHANNEL)
    # direct all error printing to stdout
    error_printer = printerTo('stdout', ERROR_CHANNEL)

    exact_image = Image('my_image.hv')
    image = exact_image.get_empty_copy()
    image.fill(1.0)

    acq_templ = AcquisitionData('my_raw_data.hs')

    # create matrix to be used by the acquisition model
    matrix = RayTracingMatrix()
    matrix.set_num_tangential_LORs(2)

    # create acquisition model
    #am = stir.AcquisitionModelUsingMatrix()
    am = AcquisitionModelUsingMatrix()
    am.set_matrix(matrix)
    am.set_up(acq_templ, exact_image)

    print('projecting image...')
    ad = am.forward(exact_image)

    # define a filter
    filter = CylindricFilter()

    # define the objective function
    obj_fun = PoissonLogLh_LinModMean_AcqModData()
    #obj_fun.set_max_segment_num_to_process(3)
    obj_fun.set_acquisition_model(am)
    obj_fun.set_acquisition_data(ad)

    num_subiterations = 2

    # create OSMAPOSL reconstructor
    recon = OSMAPOSLReconstruction()
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

    add = AcquisitionData(acq_templ)
    bkg = AcquisitionData(acq_templ)
    nrm = AcquisitionData(acq_templ)
    add.fill(args.additive)
    bkg.fill(args.additive)
    nrm.fill(args.normalisation)

    print('\n--- testing normalisation only...')
    am.set_normalisation(nrm)
    am.set_up(acq_templ, exact_image)
    print('projecting image...')
    new_ad = am.forward(exact_image)
    obj_fun.set_acquisition_model(am)
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
    obj_fun.set_acquisition_model(am)
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

# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
except error as err:
    # display error information
    print('STIR exception occured: %s' % err.value)
