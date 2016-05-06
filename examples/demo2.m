if ~libisloaded('mgadgetron')
    loadlibrary('mgadgetron')
end
if ~libisloaded('mutilities')
    loadlibrary('mutilities')
end

try
    input_data = gadgetron.MR_Acquisitions('testdata.h5');
    
    gadget1 = gadgetron.Gadget('RemoveROOversamplingGadget');
	gadget2 = gadgetron.Gadget('AcquisitionAccumulateTriggerGadget');
	gadget3 = gadgetron.Gadget('BucketToBufferGadget');
	gadget4 = gadgetron.Gadget('SimpleReconGadget');
	gadget5 = gadgetron.Gadget('ImageArraySplitGadget');
	gadget6 = gadgetron.Gadget('ExtractGadget');
	
    recon = gadgetron.ImagesReconstructor();

    recon.add_gadget('g1', gadget1);
	recon.add_gadget('g2', gadget2);
	recon.add_gadget('g3', gadget3);
	recon.add_gadget('g4', gadget4);
	recon.add_gadget('g5', gadget5);
    
    recon.set_input(input_data)
    recon.process()
    complex_images = recon.get_output();
    
    proc = gadgetron.ImagesProcessor();
    proc.add_gadget('g6', gadget6);
    images = proc.process(complex_images);

    for i = 1 : images.number()
        data = images.image_as_array(i);
        figure(1000000 + i)
        data = data/max(max(max(data)));
        imshow(data(:,:,1));
    end

    images.write('output2.h5', datestr(datetime))

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
