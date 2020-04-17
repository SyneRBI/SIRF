% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2018 - 2019 University College London
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
set_up_PET();
g.SIRF_PATH = getenv('SIRF_PATH');

try_stirtonifti(g);

function try_stirtonifti(g)
    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Starting STIR to Nifti test...')
    disp('%------------------------------------------------------------------------ %')
    
    set_up_PET();

    % Input filenames
    nifti_filename = fullfile(g.SIRF_PATH, '/data/examples/Registration/test2.nii.gz');

    % Load the image as a NiftiImageData3D
    image_nifti = sirf.Reg.NiftiImageData3D(nifti_filename);

    % Read as STIRImageData, convert to NiftiImageData3D and save to file
    image_stir = sirf.STIR.ImageData(nifti_filename);
    image_nifti_from_stir = sirf.Reg.NiftiImageData3D(image_stir);
    image_nifti_from_stir.write('results/stir_to_nifti.nii',image_nifti.get_original_datatype());

    % Compare the two
    assert(image_nifti == image_nifti_from_stir, 'Conversion from STIR to Nifti failed.');

    % Resample and then check that voxel values match
    resample = sirf.Reg.NiftyResample();
    resample.set_floating_image(image_stir);
    resample.set_reference_image(image_nifti);
    resample.set_interpolation_type_to_nearest_neighbour();
    resample.process();

    % as_array() of both original images should match
    assert(all(all(all(image_nifti.as_array() == resample.get_output().as_array()))), 'as_array() of sirf.Reg.NiftiImageData and resampled sirf.STIR.ImageData are different.')

    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Finished STIR to Nifti test.')
    disp('%------------------------------------------------------------------------ %')
end