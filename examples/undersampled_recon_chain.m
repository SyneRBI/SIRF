% Lower-level demo, 3-chain GRAPPA reconstruction of undersampled data.

if ~libisloaded('mutilities')
    fprintf('loading mutilities library...\n')
    [notfound, warnings] = loadlibrary('mutilities');
end
if ~libisloaded('mgadgetron')
    fprintf('loading mgadgetron library...\n')
    [notfound, warnings] = loadlibrary('mgadgetron');
end

try
    % define raw data source
    file = input('raw data file: ', 's');
    input_data = gadgetron.MR_Acquisitions(file);

    prep_gadgets = [{'NoiseAdjustGadget'} {'AsymmetricEchoAdjustROGadget'} ...
         {'RemoveROOversamplingGadget'}];
    preprocessed_data = input_data.process(prep_gadgets);

    % define reconstructor
    gadgets = [...
        {'AcquisitionAccumulateTriggerGadget'}, ...
        {'BucketToBufferGadget'}, ...
        {'GenericReconCartesianReferencePrepGadget'}, ...
        {'GenericReconCartesianGrappaGadget'}, ...
        {'GenericReconFieldOfViewAdjustmentGadget'}, ...
        {'GenericReconImageArrayScalingGadget'}, ...
        {'ImageArraySplitGadget'} ...
        ];
    recon = gadgetron.ImagesReconstructor(gadgets);

    % perform reconstruction
    recon.set_input(preprocessed_data)
    fprintf('reconstructing images...\n')
    recon.process()
    % get reconstructed complex images and G-factors
    complex_output = recon.get_output();
    
    % extract real images
    fprintf('processing images...\n')
    output = complex_output.real();

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
