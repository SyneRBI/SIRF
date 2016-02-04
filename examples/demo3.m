if ~libisloaded('mutilities')
    fprintf('loading mutilities library...\n')
    [notfound, warnings] = loadlibrary('mutilities');
end
if ~libisloaded('mgadgetron')
    fprintf('loading mgadgetron library...\n')
    [notfound, warnings] = loadlibrary('mgadgetron');
end

try
    % define gadgets
    gadget1 = gadgets.RemoveROOversamplingGadget();
	gadget2 = gadgets.SimpleReconGadgetSet();
	gadget3 = gadgets.ExtractGadget();
    
    % set gadgets parameters
    gadget2.set_property('trigger_dimension', 'repetition')
    gadget2.set_property('split_slices', 'true')
    
    % create reconstruction object
    recon = gadgetron.MRIReconstruction();

    % build gadgets chain
    recon.add_gadget('g1', gadget1);
	recon.add_gadget('g2', gadget2);
	recon.add_gadget('g3', gadget3);
    
    % acquisitions will be read from this HDF file
    input_data = gadgetron.ISMRMRDataset('testdata.h5');
    
    % connect to input data
    recon.set_input(input_data)
    % perform reconstruction
    recon.process()
    % get reconstructed images
    images = recon.get_output();
    
    % plot reconstructed images
    data = images.image_as_array(0);
    figure(1000000)
    data = data/max(max(max(data)));
    imshow(data(:,:,1));

    % write images to a new group in 'output3.h5'
    % named after the current date and time
    images.write('output3.h5', datestr(datetime))

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
