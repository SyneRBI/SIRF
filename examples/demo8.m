if ~libisloaded('mgadgetron')
    loadlibrary('mgadgetron')
end
if ~libisloaded('mutilities')
    loadlibrary('mutilities')
end

try
    % acquisitions will be read from this HDF file
    input_data = gadgetron.ISMRMRDAcquisitions('testdata.h5');
    fprintf('%d acquisitions found\n', input_data.number())
    
    % pre-process acquisition data
    fprintf('processing acquisitions...\n')
    interim_data = gadgetron.MR_remove_x_oversampling(input_data);
	
    % perform reconstruction
    recon = gadgetron.SimpleReconstructor();
    recon.set_input(interim_data)
    fprintf('reconstructing...\n')
    recon.process()
    interim_images = recon.get_output();
    
    % post-process reconstructed images
    fprintf('processing images...\n')
    images = gadgetron.MR_extract_real_images(interim_images);

    % create acquisition model based on the acquisition parameters
    % stored in input_data and image parameters stored in interim_images
    am = gadgetron.AcquisitionModel(input_data, interim_images);

    % use the acquisition model (forward projection) to produce acquisitions
    acqs = am.forward(interim_images);

    % compute the difference between real and modelled acquisitions
    a = -acqs.dot(input_data) / input_data.dot(input_data);
    b = 1.0;
    diff = gadgetron.AcquisitionsContainer.axpby(a, input_data, b, acqs);
    fprintf('reconstruction residual: %e\n', diff.norm()/acqs.norm())

    % apply the adjoint model (backward projection)
    imgs = am.backward(diff);

    % test that the backward projection is the adjoint of forward
    % on x = diff and y = interim_images
    fprintf('(x, F y) = %s\n', num2str(diff.dot(acqs)))
    fprintf('= (B x, y) = %s\n', num2str(imgs.dot(interim_images)))

    % test images norm and dot product
    s = imgs.norm();
    fprintf('(B x, B x) = %e = %e\n', imgs.dot(imgs), s*s)

    % test linear combination of images
    a = -1.0;
    im_diff = gadgetron.ImagesContainer.axpby(a, imgs, b, imgs);
    fprintf('0.0 = %e\n', im_diff.norm())

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
