% load STIR interface library
if ~libisloaded('mstir')
    loadlibrary('mstir')
end

try
    % direct all diagnostic printing to Matlab Command Window
    printer = stir.printerTo('stdout');
    
    %image = stir.Image('my_uniform_image_circular.hv');
    image = stir.Image();
    image.read_from_file('my_uniform_image_circular.hv')

    recon = stir.OSMAPOSLReconstruction('OSMAPOSL_test_PM_QP.par');
    
    f = stir.CylindricFilter(recon.get_inter_iteration_filter());
    f.set_strictly_less_than_radius(false)
    f.apply(image)
    f.set_strictly_less_than_radius(true)
    
    obj = stir.PoissonLogLh_LinModMean_AcqModData...
        (recon.get_objective_function());
    obj.set_sensitivity_filename('RPTsens_seg3_PM.hv')
    obj.set_input_filename('Utahscat600k_ca_seg4.hs')
    prior = obj.get_prior();
    fprintf('prior penalisation factor: %f\n', prior.get_penalisation_factor())
    prior.set_penalisation_factor(0.5)
    am = obj.get_acquisition_model();
    fprintf('tangential LORs: %d\n', am.get_matrix().get_num_tangential_LORs())
    am.get_matrix().set_num_tangential_LORs(2)

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
