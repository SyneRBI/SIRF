% Upper-level demo, illustrates the use of acquisition models and images and acquisitions algebra.
% Involves the computation of coil sensitivity maps and
% the projection from the image space into acquisition space and back.

if ~libisloaded('mgadgetron')
    loadlibrary('mgadgetron')
end
if ~libisloaded('mutilities')
    loadlibrary('mutilities')
end

try
    % acquisitions will be read from this HDF file
    input_data = gadgetron.MR_Acquisitions('testdata.h5');
    fprintf('%d acquisitions found\n', input_data.number())
    
    % pre-process acquisition data
    fprintf('processing acquisitions...\n')
    processed_data = gadgetron.MR_remove_x_oversampling(input_data);
	
    % perform reconstruction
    recon = gadgetron.MR_BasicReconstruction();
    recon.set_input(processed_data)
    fprintf('reconstructing...\n')
    recon.process()
    complex_images = recon.get_output();
    
    % post-process reconstructed images
    fprintf('processing images...\n')
    images = gadgetron.MR_extract_real_images(complex_images);

    csms = gadgetron.MR_CoilSensitivityMaps();
    fprintf('sorting acquisitions...\n')
    processed_data.sort()
    fprintf('calculating sensitivity maps...\n')
    csms.calculate(processed_data)

    % create acquisition model based on the acquisition parameters
    % stored in input_data and image parameters stored in interim_images
    am = gadgetron.MR_AcquisitionModel(processed_data, complex_images);

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
