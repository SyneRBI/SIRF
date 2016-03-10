function real_images = MR_extract_real_images(complex_images)
% removes oversampling
    handle = calllib('mgadgetron', 'mGT_newObject', ...
        'ExtractRealImagesProcessor');
    gadgetron.checkExecutionStatus('MR_extract_real_images', handle);
    real_images = gadgetron.ImagesContainer();
    real_images.handle_ = calllib('mgadgetron', 'mGT_processImages', ...
         handle, complex_images.handle_);
    gadgetron.checkExecutionStatus('MR_extract_real_images', ...
        real_images.handle_);
    calllib('mutilities', 'mDeleteObject', handle)
end