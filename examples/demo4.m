if ~libisloaded('mutilities')
    fprintf('loading mutilities library...\n')
    [notfound, warnings] = loadlibrary('mutilities');
end
if ~libisloaded('mgadgetron')
    fprintf('loading mgadgetron library...\n')
    [notfound, warnings] = loadlibrary('mgadgetron');
end

%libfunctions('mutilities')
%libfunctions('mgadgetron')

msg = 'ok';

try
    %printer = gadgetron.Printer('stdout');
    %printer = gadgetron.printerTo('stdout', 0);
    
    gadget1 = gadgets.RemoveROOversamplingGadget();
	gadget2 = gadgets.AcquisitionAccumulateTriggerGadget();
	gadget3 = gadgets.BucketToBufferGadget();
	gadget4 = gadgets.SimpleReconGadget();
	gadget5 = gadgets.ImageArraySplitGadget();
	gadget6 = gadgets.ExtractGadget();
    
    gadget2.set_property('trigger_dimension', 'repetition')
    gadget3.set_property('split_slices', 'true')
    
    recon = gadgetron.MRIReconstruction();

    recon.add_gadget('g1', gadget1);
	recon.add_gadget('g2', gadget2);
	recon.add_gadget('g3', gadget3);
	recon.add_gadget('g4', gadget4);
	recon.add_gadget('g5', gadget5);
	recon.add_gadget('g6', gadget6);
    
    input_data = gadgetron.ISMRMRDataset('testdata.h5');
    
    recon.set_input(input_data)
    recon.process()
    images = recon.get_output();
    
    data = images.image_as_array(0);
    figure(1000000)
    data = data/max(max(max(data)));
    imshow(data(:,:,1));

    images.write('output4.h5', datestr(datetime))

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
