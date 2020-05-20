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

set_up_Reg();
set_up_MR();
file = fullfile(getenv('SIRF_PATH'), 'data', 'examples', 'MR', 'zenodo', 'SIRF_recon.h5');
image_gadgetron = sirf.Gadgetron.ImageData();
image_gadgetron.read(file, 'Gadgetron', 1);
image_nifti = sirf.Reg.NiftiImageData3D(image_gadgetron);
image_nifti_from_gadgetron = sirf.Reg.NiftiImageData3D(image_gadgetron);
assert(image_nifti == image_nifti_from_gadgetron, 'Conversion from Gadgetron to Nifti failed.')
