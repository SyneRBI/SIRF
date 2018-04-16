# Imports
import os
import pSIRFReg

# Paths
SIRF_PATH     = os.environ.get('SIRF_PATH')
examples_path = SIRF_PATH + '/data/examples/Registration';
output_path   = os.getcwd() + '/results/'

# Input filenames
reference_image_filename = examples_path + "/mouseFixed.nii.gz";
floating_image_filename  = examples_path + "/mouseMoving.nii.gz";
parameter_file_aladin    = examples_path + "/paramFiles/aladin.par";
parameter_file_f3d       = examples_path + "/paramFiles/f3d.par";
matrix                   = examples_path + "/transformation_matrix.txt";
wm_im2     				 = examples_path + "/weighted_mean/regis_recon_gate2.nii";
wm_im3     				 = examples_path + "/weighted_mean/regis_recon_gate3.nii";
wm_im4     				 = examples_path + "/weighted_mean/regis_recon_gate4.nii";

# Output filenames
aladin_warped            = output_path   + "python_aladin_warped";
f3d_warped               = output_path   + "python_f3d_warped";
TM_fwrd					 = output_path   + "python_TM_fwrd.txt";
TM_back					 = output_path   + "python_TM_back.txt";
aladin_disp_fwrd 	 	 = output_path   + "python_aladin_disp_fwrd";
aladin_disp_back    	 = output_path   + "python_aladin_disp_back";
f3d_disp_fwrd  			 = output_path   + "python_f3d_disp_fwrd";
f3d_disp_back	 		 = output_path   + "python_f3d_disp_back";

output_resample          = output_path   + "python_resample";
output_activity_corr     = output_path   + "python_activity_corr";
output_weighted_mean     = output_path   + "python_weighted_mean";

# ----------------------------------------------------------------------- #
# 							Nifty aladin
#------------------------------------------------------------------------ #
NA = pSIRFReg.NiftyAladin();
NA.set_reference_image_filename      (    reference_image_filename   );
NA.set_floating_image_filename       (     floating_image_filename   );
NA.set_parameter_file		         (      parameter_file_aladin    );
NA.update();
NA.save_warped_image                 (         aladin_warped         );
NA.save_transformation_matrix_fwrd   (             TM_fwrd           );
NA.save_transformation_matrix_back   (             TM_back           );
NA.save_displacement_field_fwrd_image( aladin_disp_fwrd, True,  True );
NA.save_displacement_field_back_image( aladin_disp_back, True,  True );

# ----------------------------------------------------------------------- #
# 							Nifty f3d
#------------------------------------------------------------------------ #
NF = pSIRFReg.NiftyF3dSym();
NF.set_reference_image_filename      (  reference_image_filename  );
NF.set_floating_image_filename       (   floating_image_filename  );
NF.set_parameter_file		         (     parameter_file_f3d     );
NF.set_reference_time_point	         (             1              );
NF.set_floating_time_point	         (             1              );
NF.update();
NF.save_warped_image                 (         f3d_warped         );
NF.save_displacement_field_fwrd_image( f3d_disp_fwrd, True,  True );
NF.save_displacement_field_fwrd_image( f3d_disp_back, True,  True );

# ----------------------------------------------------------------------- #
# 							Nifty resample
#------------------------------------------------------------------------ #
NR = pSIRFReg.NiftyResample();
NR.set_reference_image_filename       ( reference_image_filename );
NR.set_floating_image_filename        ( floating_image_filename  );
NR.add_transformation_matrix_filename (          matrix          );
NR.set_interpolation_type_to_cubic_spline();
NR.update();
NR.save_resampled_image               (       output_resample    );

# ----------------------------------------------------------------------- #
# 							Activity correction
#------------------------------------------------------------------------ #
AC = pSIRFReg.ActivityCorrect();
AC.set_initial_activity    (        267000000         );
AC.set_half_life           (          6586.2          );
AC.set_input_image_filename( reference_image_filename );
AC.set_start               (            0             );
AC.set_stop                (           10             );
AC.update();
AC.save_output             (    output_activity_corr  );

# ----------------------------------------------------------------------- #
# 							Weighted mean
#------------------------------------------------------------------------ #
WM = pSIRFReg.ImageWeightedMean();
WM.add_image         (     wm_im2, 0.2    );
WM.add_image         (     wm_im3, 0.2    );
WM.add_image         (     wm_im4, 0.2    );
WM.update();
WM.save_image_to_file(output_weighted_mean);