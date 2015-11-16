% load STIR interface library
if ~libisloaded('mstir')
    loadlibrary('mstir')
end

try
    % direct all diagnostic printing to Matlab Command Window
    printer = stir.printerTo('stdout');
    
    image = stir.Image();
    image_size = [60 60 31];
    voxel_size = [4.44114 4.44114 3.375];
    image.initialise(image_size, voxel_size)    
    image.fill(1.0)

    filter = stir.CylindricFilter();
    filter.set_strictly_less_than_radius(false)
    filter.apply(image)
    filter.set_strictly_less_than_radius(true)

    data = image.density();
    figure(1000000)
    data = data/max(max(max(data)));
    imshow(data(:,:,1));

    matrix = stir.RayTracingMatrix();
    matrix.set_num_tangential_LORs(2)
    
    amd = stir.AcquisitionModelData();
    amd.read_from_file('Utahscat600k_ca_seg4.hs')

    am = stir.AcquisitionModelUsingMatrix();
    am.set_matrix(matrix);
    
    prior = stir.QuadraticPrior();
    prior.set_penalisation_factor(0.5)

    obj_fun = stir.PoissonLogLh_LinModMean_AcqModData();
    obj_fun.set_sensitivity_filename('RPTsens_seg3_PM.hv')
    obj_fun.set_use_subset_sensitivities(false)
    obj_fun.set_zero_seg0_end_planes(true)
    obj_fun.set_max_segment_num_to_process(3)
    obj_fun.set_acquisition_model(am)
    obj_fun.set_acquisition_model_data(amd)
    obj_fun.set_prior(prior)
    
    recon = stir.OSMAPOSLReconstruction();
    
    recon.set_objective_function(obj_fun)
    recon.set_MAP_model('multiplicative')
    recon.set_num_subsets(12)
    recon.set_num_subiterations(6)
    recon.set_save_interval(6)
    recon.set_inter_iteration_filter_interval(1)
    recon.set_inter_iteration_filter(filter)
    recon.set_output_filename_prefix('reconstructedImage')
    
    recon.set_up(image)

    start = recon.get_start_subiteration_num();
    stop = recon.get_num_subiterations();
    tic
    for iter = start : stop
        fprintf('\n--------------------- Subiteration %d\n',...
              recon.get_subiteration_num())
        recon.update(image)
        data = image.density();
        figure(1000000 + iter)
        data = data/max(max(max(data)));
        imshow(data(:,:,1));
    end
    toc

    expectedImage = stir.Image('test_image_PM_QP_6.hv');
    diff = expectedImage.diff_from(image);
    fprintf('difference from expected image: %e\n', diff)

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
