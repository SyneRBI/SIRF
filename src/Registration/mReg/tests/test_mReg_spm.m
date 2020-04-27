% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2018 - 2020 University College London
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

% Paths
SIRF_PATH     = getenv('SIRF_PATH');
examples_path = fullfile(SIRF_PATH, '/data/examples/Registration');
output_prefix = fullfile(pwd, '/results/');

% Input filenames
g.ref_aladin_filename = fullfile(examples_path, '/test.nii.gz');
g.spm_working_folder  = fullfile(output_prefix, 'spm_working_folder');
g.spm_working_folder2 = fullfile(output_prefix, 'spm_working_folder2');
g.ref_aladin          = sirf.Reg.NiftiImageData3D( g.ref_aladin_filename );
g.spm_to_register_ref   = fullfile(output_prefix, 'spm_to_register_ref.nii');
g.spm_to_register_flo   = fullfile(output_prefix, 'spm_to_register_flo.nii');

% You can change these when debugging
try_spm = true;

if try_spm
    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Starting SPM test...')
    disp('%------------------------------------------------------------------------ %')

    % Resample an image with NiftyResample. Register SPM, check the result

    % TM
    translations = [5,  4, -5];
    euler_angles = [5, -2, -3];

    tm = sirf.Reg.AffineTransformation(translations, euler_angles);

    niftyreg_resampler = sirf.Reg.NiftyResample();
    niftyreg_resampler.set_padding_value(0.);
    niftyreg_resampler.set_reference_image(g.ref_aladin);
    niftyreg_resampler.set_floating_image(g.ref_aladin);
    niftyreg_resampler.add_transformation(tm);
    niftyreg_resampler.set_interpolation_type_to_linear();
    floating = niftyreg_resampler.forward(g.ref_aladin);

    % Register with SPM
    spm_reg = sirf.Reg.SPMRegistration();
    spm_reg.set_reference_image(g.ref_aladin);
    spm_reg.add_floating_image(floating);
    spm_reg.add_floating_image(floating);
    spm_reg.set_working_folder(g.spm_working_folder);
    spm_reg.set_working_folder_file_overwrite(true);
    spm_reg.set_delete_temp_files(false);
    spm_reg.process();
    spm_tm = spm_reg.get_transformation_matrix_forward(2);
    spm_inv_tm = spm_tm.get_inverse();

    % Check tm roughly equals inverse TM of the resampler
    estimated_euler_angles = spm_inv_tm.get_Euler_angles();
    estimated_translations_arr = spm_inv_tm.as_array();
    estimated_translations = estimated_translations_arr(1:3,4)';

    input_euler_angles = tm.get_Euler_angles();
    input_translations = translations;

    diff_euler_angles = 100. * (input_euler_angles - estimated_euler_angles) ./ input_euler_angles;
    diff_translations = 100. * (input_translations - estimated_translations) ./ input_translations;

    disp(['Input Euler angles:              ' num2str(input_euler_angles(:)')])
    disp(['Estimated Euler angles:          ' num2str(estimated_euler_angles(:)')])
    disp(['Percentage diff in Euler angles: ' num2str(diff_euler_angles(:)')])
    disp(['Input translations:              ' num2str(input_translations(:)')])
    disp(['Estimated translations:          ' num2str(estimated_translations(:)')])
    disp(['Percentage diff in translations: ' num2str(diff_translations(:)')])
    
    % Check differences are less than 1%
    assert(all(diff_euler_angles <= 1), 'SPM registration failed (angles).')
    assert(all(diff_translations <= 1), 'SPM registration failed (translations).')
    assert(spm_reg.get_output(2) == g.ref_aladin, 'SPM registration failed (image difference).')

    g.ref_aladin.write(g.spm_to_register_ref);
    floating.write(g.spm_to_register_flo);

    % Try to register via filename
    spm_reg2 = sirf.Reg.SPMRegistration();
    spm_reg2.set_reference_image_filename(g.spm_to_register_ref);
    spm_reg2.add_floating_image_filename(g.spm_to_register_flo);
    spm_reg2.add_floating_image_filename(g.spm_to_register_flo);
    spm_reg2.set_working_folder(g.spm_working_folder2);
    spm_reg2.set_working_folder_file_overwrite(true);
    spm_reg2.set_delete_temp_files(false);
    spm_reg2.process();

    for i=1:2
        spm_reg2.get_output(i).write([output_prefix 'spm_out_' num2str(i)]);
        spm_reg2.get_displacement_field_forward(i).write([output_prefix 'spm_disp_fwd_' num2str(i)]);
        spm_reg2.get_displacement_field_inverse(i).write([output_prefix 'spm_disp_inv_' num2str(i)]);
        spm_reg2.get_deformation_field_forward(i).write([output_prefix 'spm_def_fwd_' num2str(i)]);
        spm_reg2.get_deformation_field_inverse(i).write([output_prefix 'spm_def_inv_' num2str(i)]);
        spm_reg2.get_transformation_matrix_forward(i).write([output_prefix 'spm_tm_fwd_' num2str(i)]);
        spm_reg2.get_transformation_matrix_inverse(i).write([output_prefix 'spm_tm_inv_' num2str(i)]);
    end

    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Finished SPM test.')
    disp('%------------------------------------------------------------------------ %')
end