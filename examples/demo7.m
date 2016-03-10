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
