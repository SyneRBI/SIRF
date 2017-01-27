% Upper-level demo, illustrates the use of acquisition models and images and acquisitions algebra.
% Involves the computation of coil sensitivity maps and
% the projection from the image space into acquisition space and back.

set_up_mr

try
    % acquisitions will be read from this HDF file
    file = input('raw data file: ', 's');
    input_data = AcquisitionData(file);
    fprintf('%d acquisitions found\n', input_data.number())
    
    % pre-process acquisition data
    fprintf('processing acquisitions...\n')
    processed_data = preprocess_acquisitions(input_data);
	
    % perform reconstruction
    recon = SimpleReconstruction();
    recon.set_input(processed_data)
    fprintf('reconstructing...\n')
    recon.process()
    complex_images = recon.get_output();
    
    csms = CoilSensitivityMaps();
    fprintf('sorting acquisitions...\n')
    processed_data.sort()
    fprintf('calculating sensitivity maps...\n')
    csms.calculate(processed_data)

    % create acquisition model based on the acquisition parameters
    % stored in input_data and image parameters stored in interim_images
    am = AcquisitionModel(processed_data, complex_images);

    am.set_coil_sensitivity_maps(csms)

    % use the acquisition model (forward projection) to produce acquisitions
    acqs = am.forward(complex_images);

    % compute the difference between real and modelled acquisitions
    diff = acqs - processed_data;
    fprintf('reconstruction residual: %e\n', diff.norm()/acqs.norm())

    % apply the adjoint model (backward projection)
    imgs = am.backward(processed_data);

    % test that the backward projection is the adjoint of forward
    % on x = diff and y = complex_images
    fprintf('(x, F y) = %s\n', num2str(processed_data * acqs))
    fprintf('= (B x, y) = %s\n', num2str(imgs * complex_images))

    % test images norm and dot product
    s = imgs.norm();
    fprintf('(B x, B x) = %s = %e\n', num2str(imgs * imgs), s*s)

    % test linear combination of images
    im_diff = imgs - complex_images;
    fprintf('0.0 = %e\n', im_diff.norm()/complex_images.norm())

    complex_images.show()

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
