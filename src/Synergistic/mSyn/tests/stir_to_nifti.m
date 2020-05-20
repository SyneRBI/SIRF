% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2020 Rutherford Appleton Laboratory STFC.
% Copyright 2020 University College London.
% 
% This is software developed for the Collaborative Computational
% Project in Positron Emission Tomography and Magnetic Resonance imaging
% (http://www.ccppetmr.ac.uk/).
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% http://www.apache.org/licenses/LICENSE-2.0
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

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
