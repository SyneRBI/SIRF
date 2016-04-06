function output_data = MR_remove_x_oversampling(input_data)
% removes oversampling
    handle = calllib('mgadgetron', 'mGT_newObject', ...
        'RemoveOversamplingProcessor');
    gadgetron.checkExecutionStatus('MR_remove_x_oversampling', handle);
    output_data = gadgetron.AcquisitionsContainer();
    output_data.handle_ = calllib('mgadgetron', 'mGT_processAcquisitions', ...
         handle, input_data.handle_);
    gadgetron.checkExecutionStatus('MR_remove_x_oversampling', ...
        output_data.handle_);
    calllib('mutilities', 'mDeleteObject', handle)
end