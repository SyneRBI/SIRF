if ~libisloaded('mgadgetron')
    loadlibrary('mgadgetron')
end
if ~libisloaded('mutilities')
    loadlibrary('mutilities')
end

try
    % acquisitions will be read from this HDF file
    input_data = gadgetron.ISMRMRDataset('testdata.h5');
    
    % define gadgets
	gadget1 = gadgets.RemoveROOversamplingGadget();
	gadget2 = gadgets.SimpleReconstructionGadget();
	gadget3 = gadgets.ExtractGadget();
	
    % set gadgets parameters
    gadget2.set_property('trigger_dimension', 'repetition')
    gadget2.set_property('split_slices', 'true')
    
    % create reconstruction object
    recon = gadgetron.MRIReconstruction();

    % build reconstruction chain
    recon.add_gadget('g1', gadget1);
	recon.add_gadget('g2', gadget2);
    
    % connect to input data
    recon.set_input(input_data)
    % perform reconstruction
    recon.process()
    % get reconstructed images
    images = recon.get_output();
    
    % build image post-processing chain
    proc = gadgetron.ImagesProcessor();
    proc.add_gadget('g3', gadget3);

    % post-process reconstructed images
    im2 = proc.process(images);

    % plot obtained images
    data = im2.image_as_array(0);
    figure(1000000)
    data = data/max(max(max(data)));
    imshow(data(:,:,1));

    % write images to a new group in 'output4.h5'
    % named after the current date and time
    im2.write('output4.h5', datestr(datetime))

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
