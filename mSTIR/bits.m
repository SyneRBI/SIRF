    data = image.density();
    data = data/max(max(max(data)));
    figure(1000000)
    imshow(data(:,:,1));

%opengl software

%     scale = 1.0/max(max(max(data)));
%     stir.show(data, scale, 10)
%     %drawnow

%         scale = 1.0/max(max(max(data)));
%         stir.show(data, scale, 10)

%    drawnow

    %image0.initialise(60, 60, 31, 4.44114, 4.44114, 3.375)    
    image0.fill(1.0)
    image1 = image0.clone();
    image = image0.get_empty_copy();

%    obj_fun.set_input_filename('Utahscat600k_ca_seg4.hs')

%    obj_fun.set_up()

%    recon.set_start_subset_num(0)

%    recon.set_start_subiteration_num(1)

    fprintf('subiteration range: %d to %d\n', start, stop)

    recon.set_subiteration_num(start)

    % direct all diagnostic printing to Matlab Command Window
    printer = stir.printerTo('stdout');
    
    expectedImage = stir.Image('test_image_PM_QP_6.hv');
    diff = expectedImage.diff_from(image);
    fprintf('difference from expected image: %e\n', diff)

%    image.read_from_file('my_uniform_image_circular.hv')

%     obj.set_sensitivity_filename('RPTsens_seg3_PM.hv')
%     obj.set_input_filename('Utahscat600k_ca_seg4.hs')

    f.set_strictly_less_than_radius(false)
    f.apply(image)

%     obj = stir.PoissonLogLh_LinModMean_AcqModData...
%         (recon.get_objective_function());

%    am = obj.get_acquisition_model();

%    image.read_from_file('my_image0.hv')

%     image_size = [60 60 31];
%     voxel_size = [4.44114 4.44114 3.375];

%    obj_fun.set_sensitivity_filename('RPTsens_seg3_PM.hv')
%    obj_fun.set_use_subset_sensitivities(false)

%     start = recon.get_start_subiteration_num();
%     stop = recon.get_num_subiterations();
%     tic
%    for iter = start : stop

%        data = data/max(max(max(data)));
