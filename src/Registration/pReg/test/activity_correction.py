import os

SIRF_PATH = os.environ.get('SIRF_PATH')
examples_path = SIRF_PATH + '/data/examples/Registration';

import pSIRFReg

f_input  = examples_path + '/test.nii.gz';
f_output = os.getcwd()   + '/activity_correction_PYTHON.nii';

act_corr = pSIRFReg.ActivityCorrect();
act_corr.set_initial_activity(267000000);
act_corr.set_half_life(6586.2);
act_corr.set_input_image_filename(f_input);
act_corr.set_start(0);
act_corr.set_stop(10);
act_corr.update();
act_corr.save_output(f_output);