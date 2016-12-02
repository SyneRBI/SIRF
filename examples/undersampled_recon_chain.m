% Lower-level demo, 3-chain GRAPPA reconstruction of undersampled data.

select_gadgetron

try
    % define raw data source
    file = input('raw data file: ', 's');
    input_data = AcquisitionData(file);

    prep_gadgets = [{'NoiseAdjustGadget'} {'AsymmetricEchoAdjustROGadget'} ...
         {'RemoveROOversamplingGadget'}];
    preprocessed_data = input_data.process(prep_gadgets);

    % define reconstructor
    gadgets = [...
        {'AcquisitionAccumulateTriggerGadget'}, ...
        {'BucketToBufferGadget'}, ...
        {'GenericReconCartesianReferencePrepGadget'}, ...
        {'GRAPPA:GenericReconCartesianGrappaGadget'}, ...
        {'GenericReconFieldOfViewAdjustmentGadget'}, ...
        {'GenericReconImageArrayScalingGadget'}, ...
        {'ImageArraySplitGadget'} ...
        ];
    recon = ImagesReconstructor(gadgets);
    % switch off the computation of G-factors
    recon.set_gadget_property...
        ('GRAPPA', 'send_out_gfactor', false)

    % perform reconstruction
    recon.set_input(preprocessed_data)
    fprintf('reconstructing images...\n')
    recon.process()
    % get reconstructed complex images
    output = recon.get_output();
    
    % plot reconstructed images
    output.show()

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
