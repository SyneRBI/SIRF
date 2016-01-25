if ~libisloaded('mgadgetron')
    loadlibrary('mgadgetron')
end

try
    input_data = gadgetron.ISMRMRDataset('testdata.h5');
    
	reader = gadgets.GadgetIsmrmrdAcquisitionMessageReader();
    writer = gadgets.MRIImageWriter();
	gadget1 = gadgets.RemoveROOversamplingGadget();
	gadget2 = gadgets.AcquisitionAccumulateTriggerGadget();
	gadget3 = gadgets.BucketToBufferGadget();
	gadget4 = gadgets.SimpleReconGadget();
	gadget5 = gadgets.ImageArraySplitGadget();
	gadget6 = gadgets.ExtractGadget();
	endgadget = gadgets.ImageFinishGadget();

    gc = gadgetron.GadgetChain();
	gc.addReader('r1', reader);
	gc.addWriter('w1', writer);
	gc.addGadget('g1', gadget1);
	gc.addGadget('g2', gadget2);
	gc.addGadget('g3', gadget3);
	gc.addGadget('g4', gadget4);
	gc.addGadget('g5', gadget5);
	gc.addGadget('end', endgadget);

    conn = gadgetron.ClientConnector();
    images = gadgetron.ImagesList();
	conn.register_images_receiver(images)
	conn.connect('localhost', '9002')
    conn.config_gadget_chain(gc)
    conn.send_parameters(input_data.get_header())
    conn.send_acquisitions(input_data)
    conn.disconnect()

   	gc2 = gadgetron.GadgetChain();
    im2 = gadgetron.ImagesList();
	rd2 = gadgets.MRIImageReader();

	gc2.addReader('r1', rd2);
	gc2.addWriter('w1', writer);
	gc2.addGadget('g1', gadget6);
	gc2.addGadget('end', endgadget);

    conn.register_images_receiver(im2)
    conn.connect('localhost', '9002')
    conn.config_gadget_chain(gc2)
    conn.send_images(images)
    conn.disconnect()

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
