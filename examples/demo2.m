if ~libisloaded('mgadgetron')
    loadlibrary('mgadgetron')
end
if ~libisloaded('mutilities')
    loadlibrary('mutilities')
end

try
    input_data = gadgetron.ISMRMRDataset('testdata.h5');
    
	gadget1 = gadgets.RemoveROOversamplingGadget();
	gadget2 = gadgets.AcquisitionAccumulateTriggerGadget();
	gadget3 = gadgets.BucketToBufferGadget();
	gadget4 = gadgets.SimpleReconGadget();
	gadget5 = gadgets.ImageArraySplitGadget();
	gadget6 = gadgets.ExtractGadget();
	
    recon = gadgetron.MRIReconstruction();

    recon.add_gadget('g1', gadget1);
	recon.add_gadget('g2', gadget2);
	recon.add_gadget('g3', gadget3);
	recon.add_gadget('g4', gadget4);
	recon.add_gadget('g5', gadget5);
    
    images = recon.process(input_data);
    
    proc = gadgetron.ImagesProcessor();
    proc.add_gadget('g6', gadget6);
    im2 = proc.process(images);

    data = im2.image_as_array(0);
    figure(1000000)
    data = data/max(max(max(data)));
    imshow(data(:,:,1));

    im2.write('output2.h5', datestr(datetime))

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
