set_up_Reg();
set_up_MR();
file = fullfile(getenv('SIRF_PATH'), 'data', 'examples', 'MR', 'zenodo', 'SIRF_recon.h5');
image_gadgetron = sirf.Gadgetron.ImageData();
image_gadgetron.read(file, 'Gadgetron', 1);
image_nifti = sirf.Reg.NiftiImageData3D(image_gadgetron);
image_nifti_from_gadgetron = sirf.Reg.NiftiImageData3D(image_stir);
assert(image_nifti == image_nifti_from_gadgetron, 'Conversion from Gadgetron to Nifti failed.')
