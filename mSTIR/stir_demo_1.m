% load STIR interface library
if ~libisloaded('mstir')
    loadlibrary('mstir')
end

try
    % create OSMAPOSL reconstructor
    recon = stir.OSMAPOSLReconstruction('OSMAPOSL_test_PM_QP2.par');
    
    % check/redefine some parameters
    f = stir.CylindricFilter(recon.get_inter_iteration_filter());
    f.set_strictly_less_than_radius(true)    
    obj = recon.get_objective_function();
    prior = obj.get_prior();
    fprintf('prior penalisation factor: %f\n', prior.get_penalisation_factor())
    prior.set_penalisation_factor(0.001)
    am = stir.PoissonLogLh_LinModMean_AcqModData(obj).get_acquisition_model();
    fprintf('tangential LORs: %d\n', am.get_matrix().get_num_tangential_LORs())
    am.get_matrix().set_num_tangential_LORs(2)

    % read an initial estimate for the reconstructed image from a file
    image = stir.Image('my_image0.hv');

    % set up the reconstructor
    recon.set_up(image)

    % run reconstruction
    tic
    recon.reconstruct(image)
    toc
 
    % compare the reconstructed image to the expected image
    data = image.density();
    expectedImage = stir.Image('my_image.hv');
    x_data = expectedImage.density();
    figure(1000000)
    data = data/max(max(max(x_data)));
    imshow(data(:,:,20));
    figure(1000001)
    x_data = x_data/max(max(max(x_data)));
    imshow(x_data(:,:,20));

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
