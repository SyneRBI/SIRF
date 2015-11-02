import stir

try:
    matrix = stir.RayTracingMatrix()
    matrix.set_num_tangential_LORs(20)

    projectors = stir.ProjectorsUsingMatrix()
    projectors.set_matrix(matrix)

##    mx = projectors.get_matrix()
##    n = projectors.get_matrix().get_num_tangential_LORs()
##    print('n = ', n)
##    mx.set_num_tangential_LORs(2)

    prior = stir.QuadraticPrior()
    prior.set_penalisation_factor(0.25)

    filter = stir.TruncateToCylindricalFOVImageProcessor()

    obj_fun =\
        stir.PoissonLogLikelihoodWithLinearModelForMeanAndProjData()
    obj_fun.set_input_filename('Utahscat600k_ca_seg4.hs0')
    obj_fun.set_sensitivity_filename('RPTsens_seg3_PM.hv0')
    obj_fun.set_use_subset_sensitivities(False)
    #obj_fun.set_recompute_sensitivity(True)
    obj_fun.set_zero_seg0_end_planes(True)
    obj_fun.set_max_segment_num_to_process(3)
    obj_fun.set_projector_pair(projectors)
    obj_fun.set_prior(prior)

    proj = obj_fun.get_projector_pair()
    mx = proj.get_matrix()
    n = proj.get_matrix().get_num_tangential_LORs()
    print('tangential_LORs:', n)
    mx.set_num_tangential_LORs(2)

    recon = stir.OSMAPOSLReconstruction()
    recon.set_MAP_model('multiplicative')
    recon.set_num_subsets(12)
    recon.set_start_subset_num(0)
    recon.set_num_subiterations(6)
    recon.set_start_subiteration_num(1)
    recon.set_save_interval(6)
    recon.set_inter_iteration_filter_interval(1)
    recon.set_inter_iteration_filter(filter)
    recon.set_output_filename_prefix('reconstructedImage')
    recon.set_objective_function(obj_fun)

    obj = recon.get_objective_function()
    obj.set_input_filename('Utahscat600k_ca_seg4.hs')
    obj.set_sensitivity_filename('RPTsens_seg3_PM.hv')

    prr = obj.get_prior()
    prr.set_penalisation_factor(0.5)

    image = stir.Image('my_uniform_image_circular.hv')
    #image.read_from_file('my_uniform_image_circular.hv')

    obj_fun.set_up()
    recon.set_up(image)
    recon.reconstruct(image)
    expectedImage = stir.Image('test_image_PM_QP_6.hv')
    diff = expectedImage.diff_from(image)
    print('difference from expected image:', diff)

except stir.error as err:
    print('STIR exception occured:\n', err.value)
