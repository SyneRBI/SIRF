% load STIR interface library
if ~libisloaded('mstir')
    loadlibrary('mstir')
end

try

%     % direct all information printing to a file
%     info_printer = stir.printerTo('stir_demo4_inf.txt', 0);
%     % direct all warning and error printing to stdout
%     warning_printer = stir.printerTo('stdout', 1);
%     error_printer = stir.printerTo('stdout', 2);
    
    printer = stir.Printer...
        ('stir_demo4_info.txt', 'stir_demo4_warn.txt', 'stir_demo4_errr.txt');

    image = stir.Image();
    image_size = [111, 111, 31];
    voxel_size = [3, 3, 3.375];
    image.initialise(image_size, voxel_size)

    shape = stir.EllipsoidalCylinder();

    shape.set_length(400);
    shape.set_radii([40, 100]);
    shape.set_origin([60, 0, 10]);
    image.add_shape(shape, 1.0)

    shape.set_radii([30, 30])
    shape.set_origin([-30, 60, 10])
    image.add_shape(shape, 1.5)

    shape.set_origin([-30, -60, 10])
    image.add_shape(shape, 0.75)

    z = int32(image_size(3)/2);

    data = image.density();
    figure(1)
    data = data/max(max(max(data)));
    imshow(data(:,:,z));

    % define the matrix to be used by the acquisition model
    matrix = stir.RayTracingMatrix();
    matrix.set_num_tangential_LORs(2)

    % define the acquisition model
    am = stir.AcquisitionModelUsingMatrix();
    am.set_matrix(matrix)

    % define a prior
    prior = stir.QuadraticPrior();
    prior.set_penalisation_factor(0.001)

    % define a filter
    filter = stir.CylindricFilter();

    % create an initial image estimate
    reconstructedImage = stir.Image();
    reconstructedImage.initialise(image_size, voxel_size)
    reconstructedImage.fill(1.0)
    % apply filter to get a cylindric initial image
    filter.apply(reconstructedImage)

    am.set_up('Utahscat600k_ca_seg4.hs', image)
    ad = am.forward(image, 'demo4data.hs');

    obj_fun = stir.PoissonLogLh_LinModMean_AcqModData();
    obj_fun.set_zero_seg0_end_planes(true)
    obj_fun.set_max_segment_num_to_process(3)
    obj_fun.set_acquisition_model(am)
    obj_fun.set_acquisition_data(ad)
    obj_fun.set_prior(prior)

    num_subiterations = 6;

    recon = stir.OSMAPOSLReconstruction();
    recon.set_objective_function(obj_fun)
    recon.set_MAP_model('multiplicative')
    recon.set_num_subsets(12)
    recon.set_num_subiterations(num_subiterations)
    recon.set_save_interval(num_subiterations)
    recon.set_inter_iteration_filter_interval(1)
    recon.set_inter_iteration_filter(filter)
    recon.set_output_filename_prefix('reconstructedImage')

    recon.set_up(reconstructedImage)

    data = reconstructedImage.density();
    figure(1000000)
    data = data/max(max(max(data)));
    imshow(data(:,:,z));

    for iter = 1 : num_subiterations
        fprintf('\n--------------------- Subiteration %d\n',...
              recon.get_subiteration_num())
        % perform an iteration
        recon.update(reconstructedImage)
        % plot the current image
        data = reconstructedImage.density();
        figure(1000000 + iter)
        imshow(data(:,:,z)/max(max(max(data))));
    end

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
