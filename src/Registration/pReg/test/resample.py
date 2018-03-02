import os

SIRF_PATH = os.environ.get('SIRF_PATH')
examples_path = SIRF_PATH + '/data/examples/Registration';

import pSIRFReg

reference = examples_path + "/test.nii.gz";
floating  = examples_path + "/test2.nii.gz";
matrix    = examples_path + "/transformation_matrix.txt";
output    = os.getcwd()   + "/resampled_image_MATLAB.nii";

resample = pSIRFReg.NiftyResample();
resample.set_reference_image_filename       ( reference );
resample.set_floating_image_filename        ( floating  );
resample.add_transformation_matrix_filename (  matrix   );
resample.set_interpolation_type_to_cubic_spline();
resample.update();

resample.save_resampled_image               (  output   );
