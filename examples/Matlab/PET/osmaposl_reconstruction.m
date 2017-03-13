% OSMAPOSL reconstruction demo with all parameters defined in the script
% and user-controlled iterations.

set_up_pet
import PET.*

try
    % information on computation progress to go to this file
    printer_info = printerTo('info.txt', 0);

    % warning and error messages to go to Matlab Command Window
    printer_warn = printerTo('stdout', 1);
    printer_errr = printerTo('stdout', 2);

    % create acquisition model
    am = AcquisitionModelUsingMatrix();
    
    % read acquisition model data
    [filename, pathname] = uigetfile('*.hs', 'Select raw data file');
    ad = AcquisitionData(fullfile(pathname, filename));

    % create initial image estimate
    image = ad.create_empty_image(1.0);

    % create objective function
    obj_fun = PoissonLogLh_LinModMean_AcqModData();
    obj_fun.set_acquisition_model(am)
    obj_fun.set_acquisition_data(ad)
    
    num_subiterations = 2;
    
    % create OSMAPOSL reconstructor
    recon = OSMAPOSLReconstruction();    
    recon.set_objective_function(obj_fun)
    recon.set_num_subsets(12)
    
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
        fprintf('\n--------------------- Subiteration %d\n', iter)
        % perform an iteration
        recon.update(image)
        % plot the current image
        data = image.as_array();
        figure(iter)
        imshow(data(:,:,20)/max(max(max(data))));
    end

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
