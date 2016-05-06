if ~libisloaded('mgadgetron')
    loadlibrary('mgadgetron')
end
if ~libisloaded('mutilities')
    loadlibrary('mutilities')
end

try
    % acquisitions will be read from this HDF file
    input_data = gadgetron.MR_Acquisitions('testdata.h5');
    
    % define gadgets
    gadget1 = gadgetron.Gadget('RemoveROOversamplingGadget');
    gadget2 = gadgetron.Gadget('SimpleReconGadgetSet');
    gadget3 = gadgetron.Gadget('ExtractGadget');
	
    % set gadgets parameters
    gadget2.set_property('trigger_dimension', 'repetition')
    gadget2.set_property('split_slices', 'true')
    
    % build reconstruction chain
    recon = gadgetron.ImagesReconstructor();
    recon.add_gadget('g1', gadget1);
	recon.add_gadget('g2', gadget2);
    % connect to input data
    recon.set_input(input_data)
    % perform reconstruction
    recon.process()
    % get reconstructed images
    interim_images = recon.get_output();
    
    % build image post-processing chain
    proc = gadgetron.ImagesProcessor();
    proc.add_gadget('g3', gadget3);

    % post-process reconstructed images
    images = proc.process(interim_images);

    % plot obtained images
    for i = 1 : images.number()
        data = images.image_as_array(i);
        figure(1000000 + i)
        data = data/max(max(max(data)));
        imshow(data(:,:,1));
    end

    % write images to a new group in 'output4.h5'
    % named after the current date and time
    images.write('output4.h5', datestr(datetime))

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
