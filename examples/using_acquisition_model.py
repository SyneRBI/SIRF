import argparse
import numpy
import pylab
import os
import sys
sys.path.append(os.environ.get('CSTIR_SRC') + '/../pSTIR')

from pStir import *

parser = argparse.ArgumentParser(description = \
'''
Forward projection demo: creates an image, forward-projects it to simulate
acquisition data and uses this data to reconstruct this image
''')
args = parser.parse_args()

def main():

    # output goes to files
    printer = Printer('info.txt', 'warn.txt', 'errr.txt')

    # create an empty image
    image = PETImage()
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
    am = PETAcquisitionModelUsingMatrix()
    am.set_matrix(matrix)

    # define a prior
    prior = QuadraticPrior()
    prior.set_penalisation_factor(0.001)

    # define a filter
    filter = CylindricFilter()

    # create an initial image estimate
    reconstructedImage = PETImage()
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
    # 'Utahscat600k_ca_seg4.hs' is used as a template
    templ = PETAcquisitionData('Utahscat600k_ca_seg4.hs')
    am.set_up(templ, image)
    ad = am.forward(image)
    # if the raw data is very large, it can be stored in a file
    # ad = am.forward(image, 'demo4data.hs')

    print('back-projecting image...')
    # backward-project the computed forward projection
    update = am.backward(ad)

    # define the objective function
    obj_fun = PoissonLogLh_LinModMean_AcqModData()
    obj_fun.set_max_segment_num_to_process(3)
    obj_fun.set_acquisition_model(am)
    obj_fun.set_acquisition_data(ad)
    obj_fun.set_prior(prior)

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
    recon.set_up(reconstructedImage)

    for iter in range(1, num_subiterations + 1):
        print('\n--------------------- Subiteration %d'\
              % recon.get_subiteration_num())
        # perform an iteration
        recon.update(reconstructedImage)
        # plot the current image
        data = reconstructedImage.as_array()
        pylab.figure(iter + 1)
        pylab.imshow(data[z,:,:])
        print('close Figure %d window to continue' % (iter + 1))
        pylab.show()

try:
    main()
except error as err:
    print('STIR exception occured: %s' % err.value)
