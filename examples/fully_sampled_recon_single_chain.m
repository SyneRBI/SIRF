% Lower-level interface demo, creates and runs a chain of gadgets.

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

try
    % define gadgets
    gadgets = [...
        {'RemoveROOversamplingGadget'}, ...
        {'AcquisitionAccumulateTriggerGadget'}, ...
        {'BucketToBufferGadget'}, ...
        {'SimpleReconGadget'}, ...
        {'ImageArraySplitGadget'}, ...
        {'ex:ExtractGadget'} ...
        ];
    
    % create reconstructor
    recon = gadgetron.ImagesReconstructor(gadgets);

    % change a property of the gadget labelled 'ex'
    recon.set_gadget_property('ex', 'extract_mask', 5);

    % define raw data source
    input_data = gadgetron.MR_Acquisitions('testdata.h5');    
    recon.set_input(input_data)
    % perform reconstruction
    recon.process()
    % get reconstructed images
    images = recon.get_output();
    
    % plot reconstructed images
    for i = 1 : images.number()
        data = images.image_as_array(i);
        figure(i)
        data = data/max(max(max(data)));
        imshow(data(:,:,1));
    end

    % write images to a new group in 'output1.h5'
    % named after the current date and time
    fprintf('appending output1.h5...\n')
    images.write('output1.h5', datestr(datetime))

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
