function [failed, ntests] = test_niftypet()
% PET test for niftypet projector

% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
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

% Select and import SIRF MATLAB MR package so that SIRF MR objects can be 
% created in this function without using the prefix 'MR.'


set_up_PET;

% data_path = sirf.Utilities.examples_data_path('PET');
% raw_data_file = fullfile(data_path, 'mMR', 'mMR_template_span11.hs');
template_acq_data = sirf.STIR.AcquisitionData('Siemens_mMR',11);

% Get image
image = get_image();

% Get AM
acq_model = sirf.STIR.AcquisitionModelUsingNiftyPET();
acq_model.set_cuda_verbosity(true);
acq_model.set_up(template_acq_data, image);

% % Generate test data
% simulated_acq_data = acq_model.forward(image);
% simulated_acq_data_w_noise = add_noise(simulated_acq_data,10);
% 
% obj_fun = sirf.STIR.make_Poisson_loglikelihood(template_acq_data);
% obj_fun.set_acquisition_model(acq_model);
% 
% recon = sirf.STIR.OSMAPOSLReconstructor();
% recon.set_objective_function(obj_fun);
% recon.set_num_subsets(1);
% recon.set_num_subiterations(1);
% recon.set_input(simulated_acq_data_w_noise);
% 
% print('setting up, please wait...')
% initial_estimate = image.get_uniform_copy();
% recon.set_up(initial_estimate);
% 
% print('reconstructing...')
% recon.set_current_estimate(initial_estimate);
% recon.process();
% reconstructed_im = recon.get_output();

end

function cyl = get_elliptical_cylinder(radius_x, radius_y, length, origin)
	cyl = sirf.STIR.EllipticCylinder();
    cyl.set_radii([radius_x,radius_y]);
    cyl.set_length(length);
    if exist('origin','var')
        cyl.set_origin(origin);
    end
end

function image = get_image()
    im_size = [127, 320, 320];
    im_spacing = [2.03125, 2.08626, 2.08626];
    image = sirf.STIR.ImageData();
    image.initialise(im_size, im_spacing);
    image.fill(0);
    
    cyl = get_elliptical_cylinder(200,100,1000);
    image.add_shape(cyl, 0.75);
    cyl = get_elliptical_cylinder(100,50,300,[20,30,10]);
    image.add_shape(cyl, 3);

    cyl = get_elliptical_cylinder(10,150,700,[-20,50,50]);
    image.add_shape(cyl, 1.5);
end

function proj_data_noisy = add_noise(proj_data,noise_factor)
    if ~exist('noise_factor','var')
        noise_factor = 1;
    end
    proj_data_arr = proj_data.as_array() / noise_factor;
    % Data should be >=0 anyway, but add abs just to be safe
    proj_data_arr = abs(proj_data_arr);
    proj_data_arr_noisy = imnoise(proj_data_arr,'gaussian');
    proj_data_noisy = proj_data.clone();
    proj_data_noisy.fill(proj_data_arr_noisy);
end