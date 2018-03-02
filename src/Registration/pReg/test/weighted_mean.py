import os

SIRF_PATH = os.environ.get('SIRF_PATH')
examples_path = SIRF_PATH + '/data/examples/Registration';

import pSIRFReg

im1    = examples_path + "/weighted_mean/regis_recon_gate1.nii";
im2    = examples_path + "/weighted_mean/regis_recon_gate2.nii";
im3    = examples_path + "/weighted_mean/regis_recon_gate3.nii";
im4    = examples_path + "/weighted_mean/regis_recon_gate4.nii";
output = os.getcwd()   + "/weighted_mean";

weighted_mean = pSIRFReg.ImageWeightedMean();
weighted_mean.add_image( im1, 0.2 );
weighted_mean.add_image( im2, 0.2 );
weighted_mean.add_image( im3, 0.2 );
weighted_mean.add_image( im4, 0.2 );
weighted_mean.update();
weighted_mean.save_image_to_file(output);
