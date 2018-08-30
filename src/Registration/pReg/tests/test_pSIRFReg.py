# Imports
import os
import sys
import pSIRFReg
import pSTIR
import numpy as np

# Paths
SIRF_PATH     = os.environ.get('SIRF_PATH')
examples_path = SIRF_PATH + '/data/examples/Registration'
output_path   = os.getcwd() + '/results/'

# Input filenames
reference_image_filename = examples_path + "/mouseFixed.nii.gz"
floating_image_filename  = examples_path + "/mouseMoving.nii.gz"
parameter_file_aladin    = examples_path + "/paramFiles/aladin.par"
parameter_file_f3d       = examples_path + "/paramFiles/f3d.par"
matrix                   = examples_path + "/transformation_matrix.txt"
stir_nifti               = examples_path + "/nifti_created_by_stir.nii"

# Output filenames
aladin_warped            = output_path   + "python_aladin_warped"
f3d_warped               = output_path   + "python_f3d_warped"
TM_fwrd					 = output_path   + "python_TM_fwrd.txt"
TM_back					 = output_path   + "python_TM_back.txt"
aladin_disp_fwrd 	 	 = output_path   + "python_aladin_disp_fwrd"
aladin_disp_back    	 = output_path   + "python_aladin_disp_back"
f3d_disp_fwrd  			 = output_path   + "python_f3d_disp_fwrd"
f3d_disp_back	 		 = output_path   + "python_f3d_disp_back"

output_resample          = output_path   + "python_resample"
output_activity_corr     = output_path   + "python_activity_corr"
output_weighted_mean     = output_path   + "python_weighted_mean"

output_stir_nifti        = output_path   + "python_stir_nifti.nii"

reference = pSIRFReg.ImageData( reference_image_filename )
floating  = pSIRFReg.ImageData(  floating_image_filename )
nifti     = pSIRFReg.ImageData(        stir_nifti        )

required_percentage_accuracy = float(1)

sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
sys.stderr.write(  '#                             Starting Nifty aladin test...                         #\n')
sys.stderr.write(  '# --------------------------------------------------------------------------------- #\n')
NA = pSIRFReg.NiftyAladinSym()
NA.set_reference_image               (       reference        )
NA.set_floating_image                (        floating        )
NA.set_parameter_file		         ( parameter_file_aladin  )
NA.update()
NA.save_warped_image                 (      aladin_warped     )
NA.save_transformation_matrix_fwrd   (         TM_fwrd        )
NA.save_transformation_matrix_back   (         TM_back        )
NA.save_displacement_field_fwrd_image( aladin_disp_fwrd, True )
NA.save_displacement_field_back_image( aladin_disp_back, True )
sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
sys.stderr.write(  '#                             Finished Nifty aladin test.                           #\n')
sys.stderr.write(  '# --------------------------------------------------------------------------------- #\n')




sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
sys.stderr.write(  '#                             Starting Nifty f3d test...                            #\n')
sys.stderr.write(  '# --------------------------------------------------------------------------------- #\n')
NF = pSIRFReg.NiftyF3dSym()
NF.set_reference_image               (      reference      )
NF.set_floating_image                (      floating       )
NF.set_parameter_file		         ( parameter_file_f3d  )
NF.set_reference_time_point	         (          1          )
NF.set_floating_time_point	         (          1          )
NF.update()
NF.save_warped_image                 (     f3d_warped      )
NF.save_displacement_field_fwrd_image( f3d_disp_fwrd, True )
NF.save_displacement_field_fwrd_image( f3d_disp_back, True )
sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
sys.stderr.write(  '#                             Finished Nifty f3d test.                              #\n')
sys.stderr.write(  '# --------------------------------------------------------------------------------- #\n')





sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
sys.stderr.write(  '#                             Starting Nifty resample test...                       #\n')
sys.stderr.write(  '# --------------------------------------------------------------------------------- #\n')
NR = pSIRFReg.NiftyResample()
NR.set_reference_image                (    reference    )
NR.set_floating_image                 (    floating     )
NR.set_transformation_matrix          (     matrix      )
NR.set_interpolation_type_to_cubic_spline()
NR.update()
NR.save_resampled_image               ( output_resample )
sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
sys.stderr.write(  '#                             Finished Nifty resample test.                         #\n')
sys.stderr.write(  '# --------------------------------------------------------------------------------- #\n')




sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
sys.stderr.write(  '#                             Starting weighted mean test...                        #\n')
sys.stderr.write(  '# --------------------------------------------------------------------------------- #\n')
WM = pSIRFReg.ImageWeightedMean()
WM.add_image( nifti, 0.2 )
WM.add_image( nifti, 0.2 )
WM.add_image( nifti, 0.2 )
WM.update()
WM.save_image_to_file(output_weighted_mean)
sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
sys.stderr.write(  '#                             Finished weighted mean test.                          #\n')
sys.stderr.write(  '# --------------------------------------------------------------------------------- #\n')




sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
sys.stderr.write(  '#                             Starting PET SIRFImageData test...                    #\n')
sys.stderr.write(  '# --------------------------------------------------------------------------------- #\n')
# Open stir image
pet_image_data = pSTIR.ImageData(stir_nifti)
image_data_from_stir = pSIRFReg.ImageData(pet_image_data)
# Compare to nifti IO (if they don't match, you'll see a message but don't throw an error for now)
image_data_from_nifti = pSIRFReg.ImageData(stir_nifti)
pSIRFReg.do_nifti_images_match(image_data_from_stir, image_data_from_nifti, required_percentage_accuracy)
# Print info
pSIRFReg.dump_nifti_info([image_data_from_stir, image_data_from_nifti, image_data_from_nifti])
# Save the one opened by stir
image_data_from_stir.save_to_file(output_stir_nifti)
# Now clone the converted and fill with 1's
cloned = pet_image_data
cloned.fill(1.)
# Fill the cloned image with data from converted
image_data_from_stir.copy_data_to(cloned)
sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
sys.stderr.write(  '#                             Finished PET SIRFImageData test.                      #\n')
sys.stderr.write(  '# --------------------------------------------------------------------------------- #\n')
