SIRF_PATH     = getenv('SIRF_PATH');
examples_path = [SIRF_PATH '/data/examples/Registration'];

input  = [examples_path '/test.nii.gz'];
output = [pwd           '/activity_correction_MATLAB.nii'];

% Run the test
act_corr = mSIRFReg.ActivityCorrect();
act_corr.set_initial_activity(267000000);
act_corr.set_half_life(6586.2);
act_corr.set_input_image_filename(input);
act_corr.set_start(0);
act_corr.set_stop(10);
act_corr.update();
act_corr.save_output(output);