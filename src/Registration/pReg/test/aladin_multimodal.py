import os

SIRF_PATH = os.environ.get('SIRF_PATH')
examples_path = SIRF_PATH + '/data/examples/Registration';

import pSIRFReg

reference_image_filename = examples_path + "/test.nii.gz";
floating_image_filename  = examples_path + "/standard.nii.gz";
parameter_file_aladin    = examples_path + "/paramFiles/aladin.par";
warped_image_filename    = os.getcwd()   + "/aladin_multimodal_python";

# Nifty aladin
NA = pSIRFReg.NiftyAladin();
NA.set_reference_image_filename(  reference_image_filename  )
NA.set_floating_image_filename (  floating_image_filename   )
NA.set_parameter_file          (   parameter_file_aladin    )
NA.update()
NA.save_warped_image           (   warped_image_filename    )