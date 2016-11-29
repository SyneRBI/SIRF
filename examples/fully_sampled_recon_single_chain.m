% Lower-level interface demo, creates and runs a chain of gadgets.

select_gadgetron

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
    recon = ImagesReconstructor(gadgets);

    % change a property of the gadget labelled 'ex'
    recon.set_gadget_property('ex', 'extract_mask', 5);

    % define raw data source
    input_data = AcquisitionData('testdata.h5');    
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

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
