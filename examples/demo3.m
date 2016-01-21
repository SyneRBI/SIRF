% load STIR interface library
if ~libisloaded('mgadgetron')
    loadlibrary('mgadgetron')
end

try
	reader = gadgets.GadgetIsmrmrdAcquisitionMessageReader();
    writer = gadgets.MRIImageWriter();
	gadget1 = gadgets.RemoveROOversamplingGadget();
	gadget2 = gadgets.AcquisitionAccumulateTriggerGadget();
	gadget3 = gadgets.BucketToBufferGadget();
	gadget4 = gadgets.SimpleReconGadget();
	gadget5 = gadgets.ImageArraySplitGadget();
	gadget6 = gadgets.ExtractGadget();
	gadget7 = gadgets.ImageFinishGadget();
	gc = gadgetron.GadgetChain();
	gc.addReader('r1', reader);
	gc.addWriter('w1', writer);
	gc.addGadget('g1', gadget1);
	gc.addGadget('g2', gadget2);
	gc.addGadget('g3', gadget3);
	gc.addGadget('g4', gadget4);
	gc.addGadget('g5', gadget5);
	gc.addGadget('g6', gadget6);
	gc.addGadget('g7', gadget7);

    input_data = gadgetron.ISMRMRDataset('testdata.h5');
    recon = gadgetron.MRReconstruction();
    recon.process(gc, input_data)
	images = recon.get_output();

    data = images.image_as_array(0);
    figure(1000000)
    data = data/max(max(max(data)));
    imshow(data(:,:,1));

    images.write('output3.h5', datestr(datetime))

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
