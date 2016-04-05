if ~libisloaded('mgadgetron')
    loadlibrary('mgadgetron')
end
if ~libisloaded('mutilities')
    loadlibrary('mutilities')
end

try
    % acquisitions will be read from this HDF file
    input_data = gadgetron.ISMRMRDAcquisitions('testdata.h5');
    
    % pre-process acquisition data
    fprintf('processing acquisitions...\n')
    processed_data = gadgetron.MR_remove_x_oversampling(input_data);
	
    % perform reconstruction
    recon = gadgetron.SimpleReconstructor();
    recon.set_input(processed_data)
    fprintf('reconstructing...\n')
    recon.process()
    complex_images = recon.get_output();
    
    csms = gadgetron.MRCoilSensitivityMaps();
    fprintf('ordering acquisitions...\n')
    input_data.order()
    fprintf('computing sensitivity maps...\n')
    csms.compute(input_data)

    % create acquisition model based on the acquisition parameters
    % stored in input_data and image parameters stored in interim_images
    am = gadgetron.AcquisitionModel(input_data, complex_images);

    am.set_coil_sensitivity_maps(csms)

    % use the acquisition model (forward projection) to produce acquisitions
    acqs = am.forward(complex_images);

    % compute the difference between real and modelled acquisitions
    a = -acqs.dot(input_data) / input_data.dot(input_data);
    b = 1.0;
    diff = gadgetron.AcquisitionsContainer.axpby(a, input_data, b, acqs);
    fprintf('reconstruction residual: %e\n', diff.norm()/acqs.norm())

    % apply the adjoint model (backward projection)
    imgs = am.backward(diff);

    % post-process reconstructed images
    fprintf('processing images...\n')
    images = gadgetron.MR_extract_real_images(complex_images);

    % plot obtained images
    for i = 1 : images.number()
        data = images.image_as_array(i);
        figure(1000000 + i)
        data = data/max(max(max(data)));
        imshow(data(:,:,1,1));
    end
    
catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
