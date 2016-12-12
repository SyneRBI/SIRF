% Lower-level interface demo, runs 3 gadget chains of different type:
% - acquisition processing chain,
% - reconstruction chain,
% - image processing chain

select_gadgetron

try
    % acquisitions will be read from this HDF file
    input_data = AcquisitionData('simulated_MR_2D_cartesian.h5');
    
    % process data using Acquisitions processing chain
    acq_proc = AcquisitionsProcessor({'RemoveROOversamplingGadget'});
    fprintf('processing acquisitions...\n')
    processed_data = acq_proc.process(input_data);
	
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
    img_proc = ImagesProcessor({'ExtractGadget'});
    fprintf('processing images...\n')
    images = img_proc.process(complex_images);

    % plot obtained images
    images.show()
    
catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
