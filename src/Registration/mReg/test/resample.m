SIRF_PATH     = getenv('SIRF_PATH');
examples_path = [SIRF_PATH '/data/examples/Registration'];

reference = [examples_path '/test.nii.gz'];
floating  = [examples_path '/test2.nii.gz'];
matrix    = [examples_path '/transformation_matrix.txt'];
output    = [pwd           '/resampled_image_MATLAB.nii'];

resample1 = mSIRFReg.NiftyResample();
resample1.set_reference_image_filename       ( reference );
resample1.set_floating_image_filename        ( floating  );
resample1.add_transformation_matrix_filename (  matrix   );
resample1.set_interpolation_type_to_cubic_spline();
resample1.update();

resample1.save_resampled_image               (  output   );
