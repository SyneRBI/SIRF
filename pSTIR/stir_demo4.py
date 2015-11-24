import numpy
import pylab
import stir

try:
    # create empty image
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

    # plot the image
    data = image.density()
    pylab.imshow(data[10,:,:])
    pylab.show()

    # create matrix to be used by the acquisition model
    matrix = stir.RayTracingMatrix()
    matrix.set_num_tangential_LORs(2)

    # create acquisition model
    am = stir.AcquisitionModelUsingMatrix()
    am.set_matrix(matrix)

    # create a prior
    prior = stir.QuadraticPrior()
    prior.set_penalisation_factor(0.001)

    # create filter
    filter = stir.CylindricFilter()

    # create initial image estimate
    reconstructedImage = stir.Image()
    reconstructedImage.initialise(image_size, voxel_size)
    reconstructedImage.fill(1.0)
    filter.apply(reconstructedImage)

    # forward-project the image to obtain 'raw data'
    amd = stir.AcquisitionModelData()
    amd.create_from_template_file('Utahscat600k_ca_seg4.hs')
    am.set_up(amd, image)
    am.forward(image, amd)

    # create objective function
    obj_fun = stir.PoissonLogLh_LinModMean_AcqModData()
    obj_fun.set_zero_seg0_end_planes(True)
    obj_fun.set_max_segment_num_to_process(3)
    obj_fun.set_acquisition_model(am)
    obj_fun.set_acquisition_model_data(amd)
    obj_fun.set_prior(prior)

    num_subiterations = 16

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
    pylab.imshow(data[10,:,:])
    pylab.show()

    for iter in range(1, num_subiterations + 1):
        print('\n--------------------- Subiteration ',\
              recon.get_subiteration_num())
        # perform an iteration
        recon.update(reconstructedImage)
        # plot the current image
        data = reconstructedImage.density()
        pylab.figure(iter + 1)
        pylab.imshow(data[10,:,:])
        pylab.show()

    # plot the reconstructed image
    data = reconstructedImage.density()
    pylab.figure(1)
    pylab.imshow(data[10,:,:])
    data = image.density()
    pylab.figure(2)
    pylab.imshow(data[10,:,:])
    pylab.show()

except stir.error as err:
    print('STIR exception occured:\n', err.value)
