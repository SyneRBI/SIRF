SIRF_PATH     = getenv('SIRF_PATH');
examples_path = [SIRF_PATH '/data/examples/Registration'];

reference_image_filename = [examples_path '/test.nii.gz'];
floating_image_filename  = [examples_path '/test2.nii.gz'];
parameter_file_aladin    = [examples_path '/paramFiles/aladin.par'];
warped_image_filename    = [pwd           '/aladin_longitudinal_matlab'];

%% Nifty aladin
NA = mSIRFReg.NiftyAladin();
NA.set_reference_image_filename(reference_image_filename)
NA.set_floating_image_filename (floating_image_filename )
NA.set_parameter_file          (  parameter_file_aladin )
NA.update();
NA.save_warped_image           (  warped_image_filename )