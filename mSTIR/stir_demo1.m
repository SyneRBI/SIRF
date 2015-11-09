% load STIR interface library
if ~libisloaded('mstir')
    loadlibrary('mstir')
end

opengl software

try
    % direct all diagnostic printing to Matlab Command Window
    printer = stir.printerTo('stdout');
    
    %image = stir.Image('my_uniform_image_circular.hv');
    voxel_dim = [60, 60, 31];
    voxel_size = [4.44114, 4.44114, 3.375];
    image = stir.Image();
    image.initialise(voxel_dim, voxel_size)
    image.fill(1.0)

    recon = stir.OSMAPOSLReconstruction('OSMAPOSL_test_PM_QP.par');
    
    f = stir.TruncateToCylindricalFOVImageProcessor...
        (recon.get_inter_iteration_filter());

    obj = stir.PoissonLogLikelihoodWithLinearModelForMeanAndProjData...
        (recon.get_objective_function());

%     f = recon.get_inter_iteration_filter();
%     stir.TruncateToCylindricalFOVImageProcessor(f)...
%         .set_strictly_less_than_radius(false)
    f.set_strictly_less_than_radius(false)
    f.apply(image)
    f.set_strictly_less_than_radius(true)
%     stir.TruncateToCylindricalFOVImageProcessor(f)...
%         .set_strictly_less_than_radius(true)
    
    data = image.density();
    scale = 1.0/max(max(max(data)));
    figure(1000000)
    stir.show(data, scale, 10)
    drawnow

    obj.set_sensitivity_filename('RPTsens_seg3_PM.hv')
    obj.set_input_filename('Utahscat600k_ca_seg4.hs')
    prior = obj.get_prior();
    fprintf('prior penalisation factor: %f\n', prior.get_penalisation_factor())
    prior.set_penalisation_factor(0.5)
    proj = obj.get_projector_pair();
    fprintf('tangential LORs: %d\n', proj.get_matrix().get_num_tangential_LORs())
    proj.get_matrix().set_num_tangential_LORs(2)

    tic
    obj.set_up()
    recon.set_up(image)
    recon.reconstruct(image)
    toc
 
    expectedImage = stir.Image('test_image_PM_QP_6.hv');
    diff = expectedImage.diff_from(image);
    fprintf('difference from expected image: %e\n', diff)

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
