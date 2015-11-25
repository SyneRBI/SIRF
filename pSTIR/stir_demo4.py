import numpy
import os
import pylab
import stir

try:
    # create an empty image
    image = stir.Image()
    image_size = (111, 111, 31)
    voxel_size = (3, 3, 3.375)
    image.initialise(image_size, voxel_size)

    # create a shape
    shape = stir.EllipsoidalCylinder()

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
    data = image.density()
    pylab.imshow(data[z,:,:])
    pylab.show()

    # define the matrix to be used by the acquisition model
    matrix = stir.RayTracingMatrix()
    matrix.set_num_tangential_LORs(2)

    # define the acquisition model
    am = stir.AcquisitionModelUsingMatrix()
    am.set_matrix(matrix)

    # define a prior
    prior = stir.QuadraticPrior()
    prior.set_penalisation_factor(0.001)

    # define a filter
    filter = stir.CylindricFilter()

    # create an initial image estimate
    reconstructedImage = stir.Image()
    reconstructedImage.initialise(image_size, voxel_size)
    reconstructedImage.fill(1.0)
    # apply filter to get a cylindric initial image
    filter.apply(reconstructedImage)

    # forward-project the image to obtain 'raw data'
    # 'Utahscat600k_ca_seg4.hs' is used as a template
    am.set_up('Utahscat600k_ca_seg4.hs', image)
    ad = am.forward(image, 'demo4data.hs')
    ad = stir.AcquisitionData('demo4data.hs')
    # backward-project the computed forward projection
    update = am.backward(ad)

    # define the objective function
    obj_fun = stir.PoissonLogLh_LinModMean_AcqModData()
    obj_fun.set_zero_seg0_end_planes(True)
    obj_fun.set_max_segment_num_to_process(3)
    obj_fun.set_acquisition_model(am)
    obj_fun.set_acquisition_data(ad)
    obj_fun.set_prior(prior)

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
    recon.set_up(reconstructedImage)

    # plot the initial image
    data = reconstructedImage.density()
    pylab.figure(1)
    pylab.imshow(data[z,:,:])
    pylab.show()

    for iter in range(1, num_subiterations + 1):
        print('\n--------------------- Subiteration ',\
              recon.get_subiteration_num())
        # perform an iteration
        recon.update(reconstructedImage)
        # plot the current image
        data = reconstructedImage.density()
        pylab.figure(iter + 1)
        pylab.imshow(data[z,:,:])
        pylab.show()

    # plot the reconstructed and actual images
    data = reconstructedImage.density()
    pylab.figure(1)
    pylab.imshow(data[z,:,:])
    data = image.density()
    pylab.figure(2)
    pylab.imshow(data[z,:,:])
    pylab.show()

except stir.error as err:
    print('STIR exception occured:\n', err.value)
