import os

SIRF_PATH = os.environ.get('SIRF_PATH')
examples_path = SIRF_PATH + '/data/examples/Registration';

import pSIRFReg

reference_image_filename = examples_path + "/mouseFixed.nii.gz";
floating_image_filename  = examples_path + "/mouseMoving.nii.gz";
parameter_file_aladin    = examples_path + "/paramFiles/f3d.par";
warped_image_filename    = os.getcwd()   + "/f3d_mice_python";

# Nifty aladin
NF = pSIRFReg.NiftyF3d();
NF.set_reference_image_filename(  reference_image_filename  )
NF.set_floating_image_filename (  floating_image_filename   )
NF.set_parameter_file		   (   parameter_file_aladin    )
NF.set_reference_time_point	   (              1             )
NF.set_floating_time_point	   (              1             )
NF.update()
NF.save_warped_image           (   warped_image_filename    )