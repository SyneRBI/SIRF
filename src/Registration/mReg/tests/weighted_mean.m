SIRF_PATH     = getenv('SIRF_PATH');
examples_path = [SIRF_PATH '/data/examples/Registration'];

im1    = [examples_path '/weighted_mean/regis_recon_gate1.nii'];
im2    = [examples_path '/weighted_mean/regis_recon_gate2.nii'];
im3    = [examples_path '/weighted_mean/regis_recon_gate3.nii'];
im4    = [examples_path '/weighted_mean/regis_recon_gate4.nii'];
output = [pwd           '/weighted_mean'];

% Run the test
weighted_mean1 = mSIRFReg.ImageWeightedMean();
weighted_mean1.add_image( im1, 0.2 );
weighted_mean1.add_image( im2, 0.2 );
weighted_mean1.add_image( im3, 0.2 );
weighted_mean1.add_image( im4, 0.2 );
weighted_mean1.update();
weighted_mean1.save_image_to_file(output);