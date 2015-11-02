% load STIR interface library
if ~libisloaded('mstir')
    loadlibrary('mstir')
end

opengl software

try
    % direct all diagnostic printing to Matlab Command Window
    printer = stir.printerTo('stdout');
    
    image = stir.Image();

    image.read_from_file('my_uniform_image.hv')
    data = image.density();
    nx = size(data, 1);
    ny = size(data, 2);
    nz = size(data, 3);
    scale = 1.0/max(max(max(data)));
    figure(1000001)
    stir.show(data, scale, 10)
    drawnow

    matrix = stir.RayTracingMatrix();
    %matrix.set_num_tangential_LORs(0)
    matrix.set_num_tangential_LORs(2)

    prior = stir.QuadraticPrior();
    prior.set_penalisation_factor(0.5)

    recon = stir.OSMAPOSLReconstruction('OSMAPOSL_test_PM_MRP.par');
    recon.set_output_filename_prefix('reconstructed_image')
    recon.set_MAP_model('multiplicative')
    recon.set_num_subsets(12)
    recon.set_up(image)
    recon.reconstruct(image)
 
    data = image.density();
    scale = 1.0/max(max(max(data)));
    figure(1000002)
    stir.show(data, scale, 10)
    drawnow
 
    expectedImage = stir.Image('test_image_PM_MRP_6.hv');
    diff = expectedImage.diff_from(image);
    fprintf('difference from expected image: %f\n', diff)
    
    matrix.get_num_tangential_LORs()
    prior.get_penalisation_factor()
    recon.get_num_subsets()

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
