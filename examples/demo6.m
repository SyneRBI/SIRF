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
    gadget11 = gadgetron.Gadget('NoiseAdjustGadget');
    gadget12 = gadgetron.Gadget('AsymmetricEchoGadget');
    gadget13 = gadgetron.Gadget('RemoveROOversamplingGadget');
    gadget21 = gadgetron.Gadget('AcquisitionAccumulateTriggerGadget');
    gadget22 = gadgetron.Gadget('BucketToBufferGadget');
    gadget23 = gadgetron.Gadget('PrepRefGadget');
    gadget24 = gadgetron.Gadget('CartesianGrappaGadget');
    gadget25 = gadgetron.Gadget('FOVAdjustmentGadget');
    gadget26 = gadgetron.Gadget('ScalingGadget');
    gadget27 = gadgetron.Gadget('ImageArraySplitGadget');
    gadget31 = gadgetron.Gadget('ComplexToFloatGadget');
    gadget32 = gadgetron.Gadget('FloatToShortGadget');
    
    % define raw data source
    file = input('raw data file: ', 's');
    input_data = gadgetron.MR_Acquisitions(file);

    % define acquisitions pre-processor
    acq_proc = gadgetron.AcquisitionsProcessor();
    acq_proc.add_gadget('g1', gadget11)
    acq_proc.add_gadget('g2', gadget12)
    acq_proc.add_gadget('g3', gadget13)
    fprintf('pre-processing acquisitions...\n')
    preprocessed_data = acq_proc.process(input_data);

    % define reconstructor
    recon = gadgetron.ImagesReconstructor();
    recon.add_gadget('g1', gadget21)
    recon.add_gadget('g2', gadget22)
    recon.add_gadget('g3', gadget23)
    recon.add_gadget('g4', gadget24)
    recon.add_gadget('g5', gadget25)
    recon.add_gadget('g6', gadget26)
    recon.add_gadget('g7', gadget27)    

    % perform reconstruction
    recon.set_input(preprocessed_data)
    fprintf('reconstructing images...\n')
    recon.process()
    % get reconstructed complex images and G-factors
    complex_output = recon.get_output();
    
    % extract real images
    img_proc = gadgetron.ImagesProcessor();
    img_proc.add_gadget('g1', gadget31)
    img_proc.add_gadget('g2', gadget32)
    complex_output.conversion_to_real(1)
    fprintf('processing images...\n')
    output = img_proc.process(complex_output);

    % plot reconstructed images and G-factors
    n = output.number()/2;
    fprintf('Enter slice number to view its data\n')
    fprintf('(a value outside the range [1 : %d] will stop this loop).\n', n)
    while (true)
        i = input('slice: ');
        if i < 1 || i > n
            break
        end
        data = output.image_as_array(2*i - 1);
        gdata = output.image_as_array(2*i);
        figure(i)
        data = data/max(max(max(data)));
        imshow(data(:,:,1));
        figure(i + n)
        gdata = gdata/max(max(max(gdata)));
        imshow(gdata(:,:,1));
    end

    % write images to a new group in 'output6.h5'
    % named after the current date and time
    fprintf('appending output6.h5...\n')
    output.write('output6.h5', datestr(datetime))

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
