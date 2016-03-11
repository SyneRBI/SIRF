if ~libisloaded('mgadgetron')
    loadlibrary('mgadgetron')
end
if ~libisloaded('mutilities')
    loadlibrary('mutilities')
end

try
    % acquisitions will be read from this HDF file
    input_data = gadgetron.ISMRMRDAcquisitions('testdata.h5');
    
    % define gadgets
	gadget1 = gadgets.RemoveROOversamplingGadget();
	gadget2 = gadgets.SimpleReconstructionGadget();
	gadget3 = gadgets.ExtractGadget();
    
    % build acquisitions pre-processing chain
    acq_proc = gadgetron.AcquisitionsProcessor();
    acq_proc.add_gadget('g1', gadget1)
    fprintf('processing acquisitions...\n')
    interim_data = acq_proc.process(input_data);
	
    % build reconstruction chain
    recon = gadgetron.ImagesReconstructor();
	recon.add_gadget('g2', gadget2);
    % connect to input data
    recon.set_input(interim_data)
    % perform reconstruction
    fprintf('reconstructing...\n')
    recon.process()
    % get reconstructed images
    interim_images = recon.get_output();
    
    % build image post-processing chain
    proc_img = gadgetron.ImagesProcessor();
    proc_img.add_gadget('g3', gadget3);
    % post-process reconstructed images
    fprintf('processing images...\n')
    images = proc_img.process(interim_images);

    % create acquisition model based on the acquisition parameters
    % stored in input_data and image parameters stored in interim_images
    am = gadgetron.AcquisitionModel(input_data, interim_images);

    % use the acquisition model (forward projection) to produce acquisitions
    acqs = am.forward(interim_images);
%     acqs.norm()
%     acqs.dot(input_data)

    % compute the difference between real and modelled acquisitions
    a = -acqs.dot(input_data) / input_data.dot(input_data);
    a = real(a);
    b = 1.0;
    diff = gadgetron.AcquisitionsContainer.axpby(a, input_data, b, acqs);
    fprintf('reconstruction residual: %e\n', diff.norm()/acqs.norm())

    % apply the adjoint model (backward projection)
    imgs = am.backward(acqs);

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
