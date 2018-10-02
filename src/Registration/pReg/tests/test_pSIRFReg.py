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
parameter_file_aladin = examples_path + "/paramFiles/niftyreg_aladin.par"
parameter_file_f3d = examples_path + "/paramFiles/niftyreg_f3d.par"

# Output filenames
save_nifti_image = output_path + "save_NiftiImage"
save_nifti_image_3d = output_path + "save_NiftiImage3D"
save_nifti_image_3d_tensor_not_split = output_path + "save_NiftiImage3DTensor_not_split"
save_nifti_image_3d_tensor_split = output_path + "save_NiftiImage3DTensor_split_%s"
save_nifti_image_3d_deformation_not_split = output_path + "save_NiftiImage3DDeformation_not_split"
save_nifti_image_3d_deformation_split = output_path + "save_NiftiImage3DDeformation_split_%s"
save_nifti_image_3d_displacement_not_split = output_path + "save_NiftiImage3DDisplacement_not_split"
save_nifti_image_3d_displacement_split = output_path + "save_NiftiImage3DDisplacement_split_%s"
aladin_warped = output_path + "aladin_warped"
f3d_warped = output_path + "f3d_warped"
TM_fwrd = output_path + "TM_fwrd.txt"
TM_back = output_path + "TM_back.txt"
aladin_def_fwrd = output_path + "aladin_def_fwrd"
aladin_def_back = output_path + "aladin_def_back_%s"
aladin_disp_fwrd = output_path + "aladin_disp_fwrd"
aladin_disp_back = output_path + "aladin_disp_back_%s"
f3d_def_fwrd = output_path + "f3d_disp_fwrd"
f3d_def_back = output_path + "f3d_disp_back_%s"
f3d_disp_fwrd = output_path + "f3d_disp_fwrd"
f3d_disp_back = output_path + "f3d_disp_back_%s"

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
    assert d == b, "NiftiImage deep_copy failed."

    # Addition
    e = d + d
    assert abs(e.get_max() - 2 * d.get_max()) < 0.0001, 'NiftiImage __add__/get_max() failed.'

    # Subtraction
    e = d - d
    assert abs(e.get_max()) < 0.0001, 'NiftiImage __sub__ failed.'

    # Sum
    assert abs(e.get_sum()) < 0.0001, 'NiftiImage get_sum() failed.'

    # Add num to image
    q = e + 1
    assert q.get_max() == e.get_max() + 1, 'NiftiImage __add__ val failed.'

    # Subtract num from image
    r = e - 1
    assert r.get_max() == e.get_max() - 1, 'NiftiImage __sub__ val failed.'

    # Multiply image by num
    s = e * 10
    assert s.get_max() == e.get_max() * 10, 'NiftiImage __mul__ val failed.'

    # Dimensions
    f = e.get_dimensions()
    assert np.array_equal(f, [3, 64, 64, 64, 1, 1, 1, 1]), 'NiftiImage get_dimensions() failed.'

    # Get as array
    arr = d.as_array()
    assert arr.max() == 100, 'NiftiImage as_array().max() failed.'
    assert arr.ndim == 3, 'NiftiImage as_array() ndims failed.'
    assert arr.shape == (64, 64, 64), 'NiftiImage as_array().shape failed.'

    # Test changing the datatypes
    nifti_types = list()
    nifti_types.append("NIFTI_TYPE_INT16")
    nifti_types.append("NIFTI_TYPE_INT32")
    nifti_types.append("NIFTI_TYPE_FLOAT32")
    nifti_types.append("NIFTI_TYPE_FLOAT64")
    nifti_types.append("NIFTI_TYPE_UINT8")
    nifti_types.append("NIFTI_TYPE_UINT16")
    nifti_types.append("NIFTI_TYPE_UINT32")
    nifti_types.append("NIFTI_TYPE_INT64")
    nifti_types.append("NIFTI_TYPE_UINT64")
    nifti_types.append("NIFTI_TYPE_FLOAT128")
    types = list()
    types.append('signed short')
    types.append('signed int')
    types.append('float')
    types.append('double')
    types.append('unsigned char')
    types.append('unsigned short')
    types.append('unsigned int')
    types.append('signed long long')
    types.append('unsigned long long')
    types.append('long double')

    for i in range(0, len(types)):
        aa = ref_aladin.deep_copy()
        aa.change_datatype(types[i])
        assert aa.get_max() == 255
        assert aa.get_datatype() == nifti_types[i]

    # Test dump methods
    q.dump_header()
    pSIRFReg.NiftiImage.dump_headers([q, s])

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
    assert d == b, "NiftiImage3D deep_copy failed."

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

    # Constructor from 3x3D
    d = pSIRFReg.NiftiImage3DTensor(ref_aladin, ref_aladin, ref_aladin)

    # Fill
    c.fill(100)

    # Get max
    assert c.get_max() == 100, 'NiftiImage3DTensor fill()/get_max() failed.'

    # Get min
    assert c.get_min() == 100, 'NiftiImage3DTensor fill()/get_min() failed.'

    # Deep copy
    d = c.deep_copy()
    assert d.handle != c.handle, 'NiftiImage3DTensor deep_copy failed.'
    assert d == c, "NiftiImage3DTensor deep_copy failed."

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

    # Constructor from 3x3D
    d = pSIRFReg.NiftiImage3DDisplacement(ref_aladin, ref_aladin, ref_aladin)

    # Fill
    c.fill(100)

    # Get max
    assert c.get_max() == 100, 'NiftiImage3DDisplacement fill()/get_max() failed.'

    # Get min
    assert c.get_min() == 100, 'NiftiImage3DDisplacement fill()/get_min() failed.'

    # Deep copy
    d = c.deep_copy()
    assert d.handle != c.handle, 'NiftiImage3DDisplacement deep_copy failed.'
    assert d == c, "NiftiImage3DDisplacement deep_copy failed."

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

    # Constructor from 3x3D
    d = pSIRFReg.NiftiImage3DDeformation(ref_aladin, ref_aladin, ref_aladin)

    # Fill
    c.fill(100)

    # Get max
    assert c.get_max() == 100, 'NiftiImage3DDeformation fill()/get_max() failed.'

    # Get min
    assert c.get_min() == 100, 'NiftiImage3DDeformation fill()/get_min() failed.'

    # Deep copy
    d = c.deep_copy()
    assert d.handle != c.handle, 'NiftiImage3DDeformation deep_copy failed.'
    assert d == c, "NiftiImage3DDeformation deep_copy failed."

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
    na.get_transformation_matrix_fwrd().save_to_file(TM_fwrd)
    na.get_transformation_matrix_back().save_to_file(TM_back)
    def_fwrd.save_to_file(aladin_def_fwrd)
    def_back.save_to_file_split_xyz_components(aladin_def_back)
    disp_fwrd.save_to_file(aladin_disp_fwrd)
    disp_back.save_to_file_split_xyz_components(aladin_disp_back)

    # Fwrd TM
    fwrd_tm = na.get_transformation_matrix_fwrd()
    sys.stderr.write('\nFwrd tm:\n%s\n\n' % fwrd_tm.as_array())

    # Back TM
    back_tm = na.get_transformation_matrix_back()
    sys.stderr.write('\nBack tm:\n%s\n\n' % back_tm.as_array())

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

    # Get initial transformation
    tm_init = pSIRFReg.Mat44(TM_fwrd)

    # default constructor
    nf = pSIRFReg.NiftyF3dSym()
    nf.set_reference_image(ref_f3d)
    nf.set_floating_image(flo_f3d)
    nf.set_parameter_file(parameter_file_f3d)
    nf.set_reference_time_point(1)
    nf.set_floating_time_point(1)
    nf.set_initial_affine_transformation(tm_init)
    nf.update()

    # Get outputs
    warped = nf.get_output()
    def_fwrd = nf.get_deformation_field_fwrd()
    def_back = nf.get_deformation_field_back()
    disp_fwrd = nf.get_displacement_field_fwrd()
    disp_back = nf.get_displacement_field_back()

    warped.save_to_file(f3d_warped)
    def_fwrd.save_to_file(f3d_def_fwrd)
    def_back.save_to_file_split_xyz_components(f3d_def_back)
    disp_fwrd.save_to_file(f3d_disp_fwrd)
    disp_back.save_to_file_split_xyz_components(f3d_disp_back)

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

    # Get transformations
    a3 = na.get_transformation_matrix_fwrd()
    b3 = na.get_displacement_field_fwrd()
    c3 = na.get_deformation_field_fwrd()

    # Get as deformations
    a_def = a3.get_as_deformation_field(ref_aladin)
    b_def = b3.get_as_deformation_field(ref_aladin)
    c_def = c3.get_as_deformation_field(ref_aladin)
    assert a_def == na.get_deformation_field_fwrd()
    assert b_def == na.get_deformation_field_fwrd()
    assert c_def == na.get_deformation_field_fwrd()

    # Compose into single deformation. Use two identity matrices and the disp field. Get as def and should be the same.
    tm_iden = pSIRFReg.Mat44.get_identity()
    trans = [tm_iden, tm_iden, c3]
    composed = pSIRFReg.NiftiImage3DDeformation.compose_single_deformation(trans, ref_aladin)
    assert composed == na.get_deformation_field_fwrd()

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

    tm_iden = pSIRFReg.Mat44.get_identity()
    tm      = na.get_transformation_matrix_fwrd()
    disp    = na.get_displacement_field_fwrd()
    deff    = na.get_deformation_field_fwrd()

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

    assert na.get_output() == nr1.get_output()

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
    # Change to float to avoid rounding errors
    ref_aladin_float = ref_aladin
    ref_aladin_float.change_datatype('float')
    im1 = ref_aladin_float.deep_copy()
    im2 = ref_aladin_float.deep_copy()
    im3 = ref_aladin_float.deep_copy()
    im4 = ref_aladin_float.deep_copy()
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
    res = ref_aladin_float.deep_copy()
    res.fill(4.5)
    assert wm1.get_output() == res

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
    assert wm2.get_output() == res

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
    pet_image_data = pSTIR.ImageData(ref_aladin_filename)
    image_data_from_stir = pSIRFReg.NiftiImage3D(pet_image_data)

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


# SIRFRegMat44
def try_sirfregmat44(na):
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting SIRFRegMat44 test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # Construct from file
    a = pSIRFReg.Mat44(TM_fwrd)

    # Multiply fwrd and inverse, should equal identity
    b = na.get_transformation_matrix_fwrd()
    c = na.get_transformation_matrix_back()
    d = b * c
    e = pSIRFReg.Mat44.get_identity()
    assert d == e, 'SIRFRegMat44::mult/comparison failed.'

    d.fill(3)
    f = d.as_array()
    assert np.all(f == 3), 'SIRFRegMat44::fill/operator[] failed.'

    assert d.get_determinant() < 1.e-7, 'SIRFRegMat44::get_determinant failed.'
    assert e.get_determinant() - 1. < 1.e-7, 'SIRFRegMat44::get_determinant failed.'

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished SIRFRegMat44 test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)


def test():
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
    try_sirfregmat44(na)


if __name__ == "__main__":
    test()
