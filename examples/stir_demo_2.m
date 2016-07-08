% OSMAPOSL reconstruction demo with all parameters defined in the script
% and user-controlled iterations

% load C++-to-C interface library
if ~libisloaded('mutilities')
    loadlibrary('mutilities')
end
% load STIR interface library
if ~libisloaded('mstir')
    loadlibrary('mstir')
end

try
    % information on computation progress to go to this file
    printer_info = stir.printerTo('stir_demo2.txt', 0);

    % warning and error messages to go to Matlab Command Window
    printer_warn = stir.printerTo('stdout', 1);
    printer_errr = stir.printerTo('stdout', 2);

    % create matrix to be used by the acquisition model
    matrix = stir.RayTracingMatrix();
    matrix.set_num_tangential_LORs(2)
    
    % create acquisition model
    am = stir.AcquisitionModelUsingMatrix();
    am.set_matrix(matrix);
    
    % read acquisition model data
    ad = stir.AcquisitionData();
    ad.read_from_file('my_forward_projection.hs')

    % create prior
    prior = stir.QuadraticPrior();
    prior.set_penalisation_factor(0.001)

    % create filter
    filter = stir.CylindricFilter();

    % create initial image estimate
    image = stir.Image();
    image_size = [111, 111, 31];
    voxel_size = [3, 3, 3.375];
    image.initialise(image_size, voxel_size)    
    image.fill(1.0)
    filter.set_strictly_less_than_radius(false)
    filter.apply(image)
    filter.set_strictly_less_than_radius(true)

    % create objective function
    obj_fun = stir.PoissonLogLh_LinModMean_AcqModData();
    obj_fun.set_zero_seg0_end_planes(true)
    obj_fun.set_max_segment_num_to_process(3)
    obj_fun.set_acquisition_model(am)
    obj_fun.set_acquisition_data(ad)
    obj_fun.set_prior(prior)
    
    num_subiterations = 6;
    
    %fprintf('ok\n')

    % create OSMAPOSL reconstructor
    recon = stir.OSMAPOSLReconstruction();    
    recon.set_objective_function(obj_fun)
    recon.set_MAP_model('multiplicative')
    recon.set_num_subsets(12)
    recon.set_num_subiterations(num_subiterations)
    recon.set_save_interval(num_subiterations)
    recon.set_inter_iteration_filter_interval(1)
    recon.set_inter_iteration_filter(filter)
    recon.set_output_filename_prefix('reconstructedImage')
    
    % set up the reconstructor
    recon.set_up(image)

    % plot the initial image
    data = image.as_array();
    figure(1000000)
    data = data/max(max(max(data)));
    imshow(data(:,:,1));

    % in order to see the reconstructed image evolution
    % take over the control of the iterative process
    % rather than allow recon.reconstruct to do all job at once
    for iter = 1 : num_subiterations
        fprintf('\n--------------------- Subiteration %d\n',...
              recon.get_subiteration_num())
        % perform an iteration
        recon.update(image)
        % plot the current image
        data = image.as_array();
        figure(1000000 + iter)
        imshow(data(:,:,20)/max(max(max(data))));
        % image can be post-processed
        filter.apply(image)
    end

    % compare the reconstructed image to the exact image
    exactImage = stir.Image('my_image.hv');
    x_data = exactImage.as_array();
    figure(1000000 + iter + 1)
    x_data = x_data/max(max(max(x_data)));
    imshow(x_data(:,:,20));

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
