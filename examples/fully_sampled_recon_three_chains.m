% Lower-level interface demo, runs 3 gadget chains of different type:
% - acquisition processing chain,
% - reconstruction chain,
% - image processing chain

select_gadgetron

try
    % acquisitions will be read from this HDF file
    input_data = AcquisitionData('testdata.h5');
    
    % process data using Acquisitions processing chain
    processed_data = input_data.process({'RemoveROOversamplingGadget'});
	
    % build reconstruction chain
    recon = ImagesReconstructor({'SimpleReconGadgetSet'});
    % connect to input data
    recon.set_input(processed_data)
    % perform reconstruction
    fprintf('reconstructing...\n')
    recon.process()
    % get reconstructed images
    complex_images = recon.get_output();

    % extract real images using Images processing chain
    images = complex_images.real();

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
