% load STIR interface library
if ~libisloaded('mstir')
    loadlibrary('mstir')
end

opengl software

try
    % direct all diagnostic printing to Matlab Command Window
    printer = stir.printerTo('stdout');
    
    image = stir.Image('my_uniform_image_circular.hv');
    %image.read_from_file('my_uniform_image_circular.hv')
    data = image.density();
    scale = 1.0/max(max(max(data)));
    figure(1000000)
    stir.show(data, scale, 10)
    drawnow

    matrix = stir.RayTracingMatrix();
    matrix.set_num_tangential_LORs(20)

    projectors = stir.ProjectorsUsingMatrix();
    projectors.set_matrix(matrix);
    
    prior = stir.QuadraticPrior();
    prior.set_penalisation_factor(0.05)

    filter = stir.TruncateToCylindricalFOVImageProcessor();

    obj_fun =...
        stir.PoissonLogLikelihoodWithLinearModelForMeanAndProjData();
    obj_fun.set_input_filename('Utahscat600k_ca_seg4.hs')
    obj_fun.set_sensitivity_filename('RPTsens_seg3_PM.hv')
    obj_fun.set_use_subset_sensitivities(false)
    obj_fun.set_zero_seg0_end_planes(true)
    obj_fun.set_max_segment_num_to_process(3)
    obj_fun.set_projector_pair(projectors)
    obj_fun.set_prior(prior)
    
    recon = stir.OSMAPOSLReconstruction();
    recon.set_MAP_model('multiplicative')
    recon.set_num_subsets(12)
    recon.set_start_subset_num(0)
    recon.set_num_subiterations(6)
    recon.set_save_interval(6)
    recon.set_start_subiteration_num(1)
    recon.set_inter_iteration_filter_interval(1)
    recon.set_inter_iteration_filter(filter)
    recon.set_objective_function(obj_fun)
    recon.set_output_filename_prefix('reconstructedImage')
    
    obj = recon.get_objective_function();
    prr = obj.get_prior();
    prr.set_penalisation_factor(0.5)
    proj = obj.get_projector_pair();
    mx = proj.get_matrix();
    lors = mx.get_num_tangential_LORs();
    fprintf('tangential LORs: %d\n', lors)
    mx.set_num_tangential_LORs(2)

    obj_fun.set_up()
    recon.set_up(image)

    start = recon.get_start_subiteration_num();
    stop = recon.get_num_subiterations();
    fprintf('subiteration range: %d to %d\n', start, stop)

%     tic
%     recon.reconstruct(image)
%     toc
    recon.set_subiteration_num(start)
    tic
    for iter = start : stop
        fprintf('\n--------------------- Subiteration %d\n',...
              recon.get_subiteration_num())
        recon.update(image)
%         data = image.density();
%         scale = 1.0/max(max(max(data)));
%         figure(1000000 + iter)
%         stir.show(data, scale, 10)
%         drawnow
    end
    toc

    expectedImage = stir.Image('test_image_PM_QP_6.hv');
    %expectedImage.read_from_file('test_image_PM_QP_6.hv')
    diff = expectedImage.diff_from(image);
    fprintf('difference from expected image: %e\n', diff)

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
