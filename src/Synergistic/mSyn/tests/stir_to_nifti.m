%function stir_to_nifti()
set_up_Reg();
set_up_PET();
file = fullfile(getenv('SIRF_PATH'), 'data', 'examples', 'MR', 'zenodo', 'dicom_as_nifti.nii');
image_stir = sirf.STIR.ImageData(file);
% TODO: STIR registries needed in msirf project
%image_stir = sirf.STIR.ImageData();
%image_stir.read(file, 'STIR', 1);
image_nifti = sirf.Reg.NiftiImageData3D(image_stir);
image_nifti_from_stir = sirf.Reg.NiftiImageData3D(image_stir);
assert(image_nifti == image_nifti_from_stir, 'Conversion from STIR to Nifti failed.')
