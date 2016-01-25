if ~libisloaded('mgadgetron')
    loadlibrary('mgadgetron')
end

msg = 'ok';

try
	gadget1 = gadgets.RemoveROOversamplingGadget();
	gadget2 = gadgets.AcquisitionAccumulateTriggerGadget();
	gadget3 = gadgets.BucketToBufferGadget();
	gadget4 = gadgets.SimpleReconGadget();
	gadget5 = gadgets.ImageArraySplitGadget();
	gadget6 = gadgets.ExtractGadget();
	gadget7 = gadgets.ImageFinishGadget();

    recon = gadgetron.MRIReconstruction();

    recon.addGadget('g1', gadget1);
	recon.addGadget('g2', gadget2);
	recon.addGadget('g3', gadget3);
	recon.addGadget('g4', gadget4);
	recon.addGadget('g5', gadget5);
	recon.addGadget('g6', gadget6);
	recon.addGadget('g7', gadget7);
    
    input_data = gadgetron.ISMRMRDataset('testdata.h5');
    
    recon.process(input_data)
    
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
