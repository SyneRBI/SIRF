% Lower-level interface demo, creates and runs a chain of gadgets.

select_gadgetron

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
    [filename, pathname] = uigetfile('*.h5', 'Select raw data file');
    input_data = AcquisitionData(fullfile(pathname, filename));
    recon.set_input(input_data)
    % perform reconstruction
    recon.process()
    % get reconstructed images
    images = recon.get_output();
    
    % plot obtained images
    images.show()

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
