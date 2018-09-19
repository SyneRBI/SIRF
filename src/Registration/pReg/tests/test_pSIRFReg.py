# Imports
import os
import sys
import time

import numpy as np
import pSIRFReg
import pSTIR
from pUtilities import *

# Paths
SIRF_PATH = os.environ.get('SIRF_PATH')
examples_path = SIRF_PATH + '/data/examples/Registration'
output_path = os.getcwd() + '/results/python_'

# Input filenames
ref_aladin_filename = examples_path + "/test.nii.gz"
flo_aladin_filename = examples_path + "/test2.nii.gz"
ref_f3d_filename = examples_path + "/mouseFixed.nii.gz"
flo_f3d_filename = examples_path + "/mouseMoving.nii.gz"
parameter_file_aladin = examples_path + "/paramFiles/aladin.par"
parameter_file_f3d = examples_path + "/paramFiles/f3d.par"
matrix = examples_path + "/transformation_matrix.txt"
stir_nifti = examples_path + "/nifti_created_by_stir.nii"

# Output filenames
save_nifti_image = output_path + "save_NiftiImage"
save_nifti_image_3d = output_path + "save_NiftiImage3D"
save_nifti_image_3d_tensor_not_split = output_path + "save_NiftiImage3DTensor_not_split"
save_nifti_image_3d_tensor_split = output_path + "save_NiftiImage3DTensor_split"
save_nifti_image_3d_deformation_not_split = output_path + "save_NiftiImage3DDeformation_not_split"
save_nifti_image_3d_deformation_split = output_path + "save_NiftiImage3DDeformation_split"
save_nifti_image_3d_displacement_not_split = output_path + "save_NiftiImage3DDisplacement_not_split"
save_nifti_image_3d_displacement_split = output_path + "save_NiftiImage3DDisplacement_split"
aladin_warped = output_path + "aladin_warped"
f3d_warped = output_path + "f3d_warped"
TM_fwrd = output_path + "TM_fwrd.txt"
TM_back = output_path + "TM_back.txt"
aladin_def_fwrd = output_path + "aladin_def_fwrd"
aladin_def_back = output_path + "aladin_def_back"
aladin_disp_fwrd = output_path + "aladin_disp_fwrd"
aladin_disp_back = output_path + "aladin_disp_back"
f3d_def_fwrd = output_path + "f3d_disp_fwrd"
f3d_def_back = output_path + "f3d_disp_back"
f3d_disp_fwrd = output_path + "f3d_disp_fwrd"
f3d_disp_back = output_path + "f3d_disp_back"

rigid_resample = output_path + "rigid_resample"
nonrigid_resample_disp = output_path + "nonrigid_resample_disp"
nonrigid_resample_def = output_path + "nonrigid_resample_def"
output_weighted_mean = output_path + "weighted_mean"
output_weighted_mean_def = output_path + "weighted_mean_def"

output_stir_nifti = output_path + "stir_nifti.nii"

ref_aladin = pSIRFReg.NiftiImage3D(ref_aladin_filename)
flo_aladin = pSIRFReg.NiftiImage3D(flo_aladin_filename)
ref_f3d = pSIRFReg.NiftiImage3D(ref_f3d_filename)
flo_f3d = pSIRFReg.NiftiImage3D(flo_f3d_filename)
nifti = pSIRFReg.NiftiImage3D(stir_nifti)

required_percentage_accuracy = float(1)


# Misc functions
def try_misc_functions():
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting misc functions test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # do nifti images match?
    assert pSIRFReg.do_nifti_images_match(ref_aladin, ref_aladin, required_percentage_accuracy), \
        "Images don\'t match, but they should."
    sys.stderr.write('\nThe following images intentionally do not match\n')
    assert not pSIRFReg.do_nifti_images_match(ref_aladin, flo_aladin, required_percentage_accuracy), \
        "Images match, but they shouldn't."

    # dump from filename
    pSIRFReg.dump_nifti_info(ref_aladin_filename)
    # dump from NiftiImage3D
    pSIRFReg.dump_nifti_info(ref_aladin)
    # dump from multiple images
    pSIRFReg.dump_nifti_info([ref_aladin, flo_aladin, nifti])
    # dump from NiftiImage3DDeformation
    deform = pSIRFReg.NiftiImage3DDeformation()
    deform.create_from_3D_image(ref_aladin)
    pSIRFReg.dump_nifti_info(deform)

    # identity matrix
    tm_iden = pSIRFReg.get_matrix()
    print (tm_iden)

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished misc functions test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)


# NiftiImage
def try_niftiimage():
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting NiftiImage test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # default constructor
    a = pSIRFReg.NiftiImage()

    # Read from file
    b = pSIRFReg.NiftiImage(ref_aladin_filename)

    # Save to file
    b.save_to_file(save_nifti_image)

    # Fill
    b.fill(100)

    # Get max
    assert b.get_max() == 100, 'NiftiImage fill()/get_max() failed.'

    # Get min
    assert b.get_min() == 100, 'NiftiImage fill()/get_min() failed.'

    # Deep copy
    d = b.deep_copy()
    assert d.handle != b.handle, 'NiftiImage deep_copy failed.'
    assert pSIRFReg.do_nifti_images_match(d, b, required_percentage_accuracy), "NiftiImage deep_copy failed."

    # Addition
    e = d + d
    assert abs(e.get_max() - 2 * d.get_max()) < 0.0001, 'NiftiImage __add__/get_max() failed.'

    # Subtraction
    e = d - d
    assert abs(e.get_max()) < 0.0001, 'NiftiImage __sub__ failed.'

    # Sum
    assert abs(e.get_sum()) < 0.0001, 'NiftiImage get_sum() failed.'

    # Dimensions
    f = e.get_dimensions()
    assert np.array_equal(f, [3, 64, 64, 64, 1, 1, 1, 1]), 'NiftiImage get_dimensions() failed.'

    # Get as array
    arr = d.as_array()
    assert arr.max() == 100, 'NiftiImage as_array().max() failed.'
    assert arr.ndim == 3, 'NiftiImage as_array() ndims failed.'
    assert arr.shape == (64, 64, 64), 'NiftiImage as_array().shape failed.'

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished NiftiImage test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

# NiftiImage3D
def try_niftiimage3d():
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting NiftiImage3D test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # default constructor
    a = pSIRFReg.NiftiImage3D()

    # Read from file
    b = pSIRFReg.NiftiImage3D(ref_aladin_filename)

    # Save to file
    b.save_to_file(save_nifti_image_3d)

    # Fill
    b.fill(100)

    # Get max
    assert b.get_max() == 100, 'NiftiImage3D fill()/get_max() failed.'

    # Get min
    assert b.get_min() == 100, 'NiftiImage3D fill()/get_min() failed.'

    # Deep copy
    d = b.deep_copy()
    assert d.handle != b.handle, 'NiftiImage3D deep_copy failed.'
    assert pSIRFReg.do_nifti_images_match(d, b, required_percentage_accuracy), "NiftiImage3D deep_copy failed."

    # Addition
    e = d + d
    assert abs(e.get_max() - 2 * d.get_max()) < 0.0001, 'NiftiImage3D __add__/get_max() failed.'

    # Subtraction
    e = d - d
    assert abs(e.get_max()) < 0.0001, 'NiftiImage3D __sub__ failed.'

    # Sum
    assert abs(e.get_sum()) < 0.0001, 'NiftiImage3D get_sum() failed.'

    # Dimensions
    f = e.get_dimensions()
    assert np.array_equal(f, [3, 64, 64, 64, 1, 1, 1, 1]), 'NiftiImage3D get_dimensions() failed.'

    # Get as array
    arr = d.as_array()
    assert arr.max() == 100, 'NiftiImage3D as_array().max() failed.'
    assert arr.ndim == 3, 'NiftiImage3D as_array() ndims failed.'
    assert arr.shape == (64, 64, 64), 'NiftiImage3D as_array().shape failed.'

    # Construct from stir pSTIR.ImageData
    stir = pSTIR.ImageData(stir_nifti)
    c = pSIRFReg.NiftiImage3D(stir)
    c.fill(100)

    # Copy data to pSTIRImageData
    stir.fill(3.)
    c.copy_data_to(stir)
    assert stir.as_array().max() == 100, 'NiftiImage3D copy_data_to stir ImageData failed.'

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished NiftiImage3D test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)


# NiftiImage3DTensor
def try_niftiimage3dtensor():
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting NiftiImage3DTensor test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # Create NiftiImage3DTensor from NiftiImage3D
    b = pSIRFReg.NiftiImage3DTensor()
    b.create_from_3D_image(ref_aladin)

    # # Save to file
    b.save_to_file(save_nifti_image_3d_tensor_not_split)
    b.save_to_file_split_xyz_components(save_nifti_image_3d_tensor_split)

    # Constructor from file
    c = pSIRFReg.NiftiImage3DTensor(save_nifti_image_3d_tensor_not_split + ".nii")

    # Fill
    c.fill(100)

    # Get max
    assert c.get_max() == 100, 'NiftiImage3DTensor fill()/get_max() failed.'

    # Get min
    assert c.get_min() == 100, 'NiftiImage3DTensor fill()/get_min() failed.'

    # Deep copy
    d = c.deep_copy()
    assert d.handle != c.handle, 'NiftiImage3DTensor deep_copy failed.'
    assert pSIRFReg.do_nifti_images_match(d, c, required_percentage_accuracy), "NiftiImage3DTensor deep_copy failed."

    # Addition
    e = d + d
    assert abs(e.get_max() - 2 * d.get_max()) < 0.0001, 'NiftiImage3DTensor __add__/get_max() failed.'

    # Subtraction
    e = d - d
    assert abs(e.get_max()) < 0.0001, 'NiftiImage3DTensor __sub__ failed.'

    # Sum
    assert abs(e.get_sum()) < 0.0001, 'NiftiImage3DTensor get_sum() failed.'

    # Dimensions
    f = e.get_dimensions()
    assert np.array_equal(f, [5, 64, 64, 64, 1, 3, 1, 1]), 'NiftiImage3DTensor get_dimensions() failed.'

    # Get as array
    arr = d.as_array()
    assert arr.max() == 100, 'NiftiImage3DTensor as_array().max() failed.'
    assert arr.ndim == 5, 'NiftiImage3DTensor as_array() ndims failed.'
    assert arr.shape == (64, 64, 64, 1, 3), 'NiftiImage3DTensor as_array().shape failed.'

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished NiftiImage3DTensor test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)


# NiftiImage3DDisplacement
def try_niftiimage3ddisplacement():
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting NiftiImage3DDisplacement test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # Create NiftiImage3DDisplacement from NiftiImage3D
    b = pSIRFReg.NiftiImage3DDisplacement()
    b.create_from_3D_image(ref_aladin)

    # Save to file
    b.save_to_file(save_nifti_image_3d_displacement_not_split)
    b.save_to_file_split_xyz_components(save_nifti_image_3d_displacement_split)

    # Constructor from file
    c = pSIRFReg.NiftiImage3DDisplacement(save_nifti_image_3d_displacement_not_split + ".nii")

    # Fill
    c.fill(100)

    # Get max
    assert c.get_max() == 100, 'NiftiImage3DDisplacement fill()/get_max() failed.'

    # Get min
    assert c.get_min() == 100, 'NiftiImage3DDisplacement fill()/get_min() failed.'

    # Deep copy
    d = c.deep_copy()
    assert d.handle != c.handle, 'NiftiImage3DDisplacement deep_copy failed.'
    assert pSIRFReg.do_nifti_images_match(d, c, required_percentage_accuracy), "NiftiImage3DDisplacement deep_copy failed."

    # Addition
    e = d + d
    assert abs(e.get_max() - 2 * d.get_max()) < 0.0001, 'NiftiImage3DDisplacement __add__/get_max() failed.'

    # Subtraction
    e = d - d
    assert abs(e.get_max()) < 0.0001, 'NiftiImage3DDisplacement __sub__ failed.'

    # Sum
    assert abs(e.get_sum()) < 0.0001, 'NiftiImage3DDisplacement get_sum() failed.'

    # Dimensions
    f = e.get_dimensions()
    assert np.array_equal(f, [5, 64, 64, 64, 1, 3, 1, 1]), 'NiftiImage3DDisplacement get_dimensions() failed.'

    # Get as array
    arr = d.as_array()
    assert arr.max() == 100, 'NiftiImage3DDisplacement as_array().max() failed.'
    assert arr.ndim == 5, 'NiftiImage3DDisplacement as_array() ndims failed.'
    assert arr.shape == (64, 64, 64, 1, 3), 'NiftiImage3DDisplacement as_array().shape failed.'

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished NiftiImage3DDisplacement test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)


# NiftiImage3DDeformation
def try_niftiimage3ddeformation():
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting NiftiImage3DDeformation test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # Create NiftiImage3DDeformation from NiftiImage3D
    b = pSIRFReg.NiftiImage3DDeformation()
    b.create_from_3D_image(ref_aladin)

    # Save to file
    b.save_to_file(save_nifti_image_3d_deformation_not_split)
    b.save_to_file_split_xyz_components(save_nifti_image_3d_deformation_split)

    # Constructor from file
    c = pSIRFReg.NiftiImage3DDeformation(save_nifti_image_3d_deformation_not_split + ".nii")

    # Fill
    c.fill(100)

    # Get max
    assert c.get_max() == 100, 'NiftiImage3DDeformation fill()/get_max() failed.'

    # Get min
    assert c.get_min() == 100, 'NiftiImage3DDeformation fill()/get_min() failed.'

    # Deep copy
    d = c.deep_copy()
    assert d.handle != c.handle, 'NiftiImage3DDeformation deep_copy failed.'
    assert pSIRFReg.do_nifti_images_match(d, c, required_percentage_accuracy), "NiftiImage3DDeformation deep_copy failed."

    # Addition
    e = d + d
    assert abs(e.get_max() - 2 * d.get_max()) < 0.0001, 'NiftiImage3DDeformation __add__/get_max() failed.'

    # Subtraction
    e = d - d
    assert abs(e.get_max()) < 0.0001, 'NiftiImage3DDeformation __sub__ failed.'

    # Sum
    assert abs(e.get_sum()) < 0.0001, 'NiftiImage3DDeformation get_sum() failed.'

    # Dimensions
    f = e.get_dimensions()
    assert np.array_equal(f, [5, 64, 64, 64, 1, 3, 1, 1]), 'NiftiImage3DDeformation get_dimensions() failed.'

    # Get as array
    arr = d.as_array()
    assert arr.max() == 100, 'NiftiImage3DDeformation as_array().max() failed.'
    assert arr.ndim == 5, 'NiftiImage3DDeformation as_array() ndims failed.'
    assert arr.shape == (64, 64, 64, 1, 3), 'NiftiImage3DDeformation as_array().shape failed.'

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished NiftiImage3DDeformation test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)


# Nifty aladin
def try_niftyaladin():
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting Nifty aladin test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # default constructor
    na = pSIRFReg.NiftyAladinSym()
    na.set_reference_image(ref_aladin)
    na.set_floating_image(flo_aladin)
    na.set_parameter_file(parameter_file_aladin)
    na.update()

    # Get outputs
    warped = na.get_output()
    def_fwrd = na.get_deformation_field_fwrd()
    def_back = na.get_deformation_field_back()
    disp_fwrd = na.get_displacement_field_fwrd()
    disp_back = na.get_displacement_field_back()

    warped.save_to_file(aladin_warped)
    na.save_transformation_matrix_fwrd(TM_fwrd)
    na.save_transformation_matrix_back(TM_back)
    def_fwrd.save_to_file(aladin_def_fwrd)
    def_back.save_to_file_split_xyz_components(aladin_def_back)
    disp_fwrd.save_to_file(aladin_disp_fwrd)
    disp_back.save_to_file_split_xyz_components(aladin_disp_back)

    # Fwrd TM
    fwrd_tm = na.get_transformation_matrix_fwrd()
    sys.stderr.write('\nFwrd tm:\n%s\n\n' % fwrd_tm)

    # Back TM
    back_tm = na.get_transformation_matrix_back()
    sys.stderr.write('\nBack tm:\n%s\n\n' % back_tm)

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished Nifty aladin test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    return na


# Nifty f3d
def try_niftyf3d():
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting Nifty f3d test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # default constructor
    nf = pSIRFReg.NiftyF3dSym()
    nf.set_reference_image(ref_f3d)
    nf.set_floating_image(flo_f3d)
    nf.set_parameter_file(parameter_file_f3d)
    nf.set_reference_time_point(1)
    nf.set_floating_time_point(1)
    nf.set_initial_affine_transformation(TM_fwrd)
    nf.update()

    # Get outputs
    warped = nf.get_output()
    def_fwrd = nf.get_deformation_field_fwrd()
    def_back = nf.get_deformation_field_back()
    disp_fwrd = nf.get_displacement_field_fwrd()
    disp_back = nf.get_displacement_field_back()

    warped.save_to_file(f3d_warped)
    def_fwrd.save_to_file_split_xyz_components(f3d_def_fwrd)
    def_back.save_to_file(f3d_def_back)
    disp_fwrd.save_to_file_split_xyz_components(f3d_disp_fwrd)
    disp_back.save_to_file(f3d_disp_back)

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished Nifty f3d test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)


# Transformation
def try_transformations(na):
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting Transformation test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # Affine
    sys.stderr.write('Testing affine...\n')
    a1 = pSIRFReg.TransformationAffine()
    a2 = pSIRFReg.TransformationAffine(TM_fwrd)
    a3 = pSIRFReg.TransformationAffine(na.get_transformation_matrix_fwrd())

    # Displacement
    sys.stderr.write('Testing displacement...\n')
    b1 = pSIRFReg.TransformationDisplacement()
    b2 = pSIRFReg.TransformationDisplacement(aladin_disp_fwrd + ".nii")
    b3 = pSIRFReg.TransformationDisplacement(na.get_displacement_field_fwrd())

    # Deformation
    sys.stderr.write('Testing deformation...\n')
    c1 = pSIRFReg.TransformationDeformation()
    c2 = pSIRFReg.TransformationDeformation(aladin_def_fwrd + ".nii")
    c3 = pSIRFReg.TransformationDeformation(na.get_deformation_field_fwrd())

    # Get as deformations
    a_def = a3.get_as_deformation_field(ref_aladin)
    b_def = b3.get_as_deformation_field(ref_aladin)
    c_def = c3.get_as_deformation_field(ref_aladin)
    assert pSIRFReg.do_nifti_images_match(a_def, na.get_deformation_field_fwrd(), required_percentage_accuracy)
    assert pSIRFReg.do_nifti_images_match(b_def, na.get_deformation_field_fwrd(), required_percentage_accuracy)
    assert pSIRFReg.do_nifti_images_match(c_def, na.get_deformation_field_fwrd(), required_percentage_accuracy)

    # Compose into single deformation. Use two identity matrices and the disp field. Get as def and should be the same.
    tm_iden = pSIRFReg.get_matrix()
    trans_aff_iden = pSIRFReg.TransformationAffine(tm_iden)
    trans = [trans_aff_iden, trans_aff_iden, c3]
    composed = pSIRFReg.compose_transformations_into_single_deformation(trans, ref_aladin)
    assert pSIRFReg.do_nifti_images_match(composed.get_as_deformation_field(ref_aladin),
                                          na.get_deformation_field_fwrd(), required_percentage_accuracy)

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished Transformation test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)


# Resample
def try_resample(na):
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting Nifty resample test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    tm_iden = pSIRFReg.get_matrix()
    tm_iden = pSIRFReg.TransformationAffine(tm_iden)
    tm      = pSIRFReg.TransformationAffine(na.get_transformation_matrix_fwrd())
    disp    = pSIRFReg.TransformationDisplacement(na.get_displacement_field_fwrd())
    deff    = pSIRFReg.TransformationDeformation(na.get_deformation_field_fwrd())

    sys.stderr.write('Testing rigid resample...\n')
    nr1 = pSIRFReg.NiftyResample()
    nr1.set_reference_image(ref_aladin)
    nr1.set_floating_image(flo_aladin)
    nr1.set_interpolation_type_to_cubic_spline()  # try different interpolations
    nr1.set_interpolation_type(3)  # try different interpolations (cubic)
    nr1.add_transformation_affine(tm_iden)
    nr1.add_transformation_affine(tm)
    nr1.update()
    nr1.get_output().save_to_file(rigid_resample)

    sys.stderr.write('Testing non-rigid displacement...\n')
    nr2 = pSIRFReg.NiftyResample()
    nr2.set_reference_image(ref_aladin)
    nr2.set_floating_image(flo_aladin)
    nr2.set_interpolation_type_to_sinc()  # try different interpolations
    nr2.set_interpolation_type_to_linear()  # try different interpolations
    nr2.add_transformation_disp(disp)
    nr2.update()
    nr2.get_output().save_to_file(nonrigid_resample_disp)

    sys.stderr.write('Testing non-rigid deformation...\n')
    nr3 = pSIRFReg.NiftyResample()
    nr3.set_reference_image(ref_aladin)
    nr3.set_floating_image(flo_aladin)
    nr3.set_interpolation_type_to_nearest_neighbour()  # try different interpolations
    nr3.add_transformation_def(deff)
    nr3.set_interpolation_type_to_linear()
    nr3.update()
    nr3.get_output().save_to_file(nonrigid_resample_def)

    assert pSIRFReg.do_nifti_images_match(na.get_output(), nr1.get_output(), required_percentage_accuracy)

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished Nifty resample test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)


# Weighted mean
def try_weighted_mean(na):
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting weighted mean test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # Do 3D
    wm1 = pSIRFReg.ImageWeightedMean()
    im1 = pSIRFReg.NiftiImage3D(stir_nifti)
    im2 = pSIRFReg.NiftiImage3D(stir_nifti)
    im3 = pSIRFReg.NiftiImage3D(stir_nifti)
    im4 = pSIRFReg.NiftiImage3D(stir_nifti)
    im1.fill(1)
    im2.fill(4)
    im3.fill(7)
    im4.fill(6)
    wm1.add_image(im1, 2)
    wm1.add_image(im2, 4)
    wm1.add_image(im3, 3)
    wm1.add_image(im4, 1)
    wm1.update()
    wm1.get_output().save_to_file(output_weighted_mean)
    # Answer should be 4.5, so compare it to that!
    res = pSIRFReg.NiftiImage3D(stir_nifti)
    res.fill(4.5)
    assert pSIRFReg.do_nifti_images_match(wm1.get_output(), res, required_percentage_accuracy)

    # Do 4D
    wm2 = pSIRFReg.ImageWeightedMean()
    im1 = na.get_deformation_field_fwrd().deep_copy()
    im2 = na.get_deformation_field_fwrd().deep_copy()
    im3 = na.get_deformation_field_fwrd().deep_copy()
    im4 = na.get_deformation_field_fwrd().deep_copy()
    im1.fill(1)
    im2.fill(4)
    im3.fill(7)
    im4.fill(6)
    wm2.add_image(im1, 2)
    wm2.add_image(im2, 4)
    wm2.add_image(im3, 3)
    wm2.add_image(im4, 1)
    wm2.update()
    wm2.get_output().save_to_file(output_weighted_mean_def)
    # Answer should be 4.5, so compare it to that!
    res = na.get_deformation_field_fwrd().deep_copy()
    res.fill(4.5)
    assert pSIRFReg.do_nifti_images_match(wm2.get_output(), res, required_percentage_accuracy)

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished weighted mean test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)


# STIR to SIRFReg
def try_stir_to_sirfreg():
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting STIR to SIRFReg test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # Open stir image
    pet_image_data = pSTIR.ImageData(stir_nifti)
    image_data_from_stir = pSIRFReg.NiftiImage3D(pet_image_data)

    # Compare to nifti IO (if they don't match, you'll see a message but don't throw an error for now)
    image_data_from_nifti = pSIRFReg.NiftiImage3D(stir_nifti)
    pSIRFReg.do_nifti_images_match(image_data_from_stir, image_data_from_nifti, required_percentage_accuracy)

    # Now fill the stir and sirfreg images with 1 and 100, respectively
    pet_image_data.fill(1.)
    image_data_from_stir.fill(100)
    assert(pet_image_data.as_array().max() != image_data_from_stir.get_max())

    # Fill the stir image with the sirfreg
    image_data_from_stir.copy_data_to(pet_image_data)
    assert(pet_image_data.as_array().max() == image_data_from_stir.get_max())

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished STIR to SIRFReg test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)


def test():
    try_misc_functions()
    try_niftiimage()
    try_niftiimage3d()
    try_niftiimage3dtensor()
    try_niftiimage3ddisplacement()
    try_niftiimage3ddeformation()
    na = try_niftyaladin()
    try_niftyf3d()
    try_transformations(na)
    try_resample(na)
    try_weighted_mean(na)
    try_stir_to_sirfreg()


if __name__ == "__main__":
    test()
