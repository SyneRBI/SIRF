function grappa_and_steepest_descent(engine)
% GRAPPA reconstruction with the steepest descent step
% to illustrate the use of Acquisition Model projections.

% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC.
% Copyright 2015 - 2017 University College London.
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
if nargin < 1
    engine = [];
end
import_str = setup_MR(engine);
eval(import_str)

% define raw data source
[filename, pathname] = uigetfile('*.h5', 'Select raw data file', mr_data_path);
acq_data = AcquisitionData(fullfile(pathname, filename));

% pre-process acquisitions
fprintf('---\n preprocessing...\n');
preprocessed_data = preprocess_acquisition_data(acq_data);
pd_norm = preprocessed_data.norm();

% perform reconstruction
recon = CartesianGRAPPAReconstructor();
recon.compute_gfactors(false);
recon.set_input(preprocessed_data);
fprintf('---\n reconstructing...\n');
recon.process();
image_data = recon.get_output();

% compute coil sensitivity maps
csms = CoilSensitivityData();
fprintf('---\n sorting acquisitions...\n')
preprocessed_data.sort()
fprintf('---\n calculating sensitivity maps...\n')
csms.calculate(preprocessed_data)

% create acquisition model based on the acquisition parameters
% stored in preprocessed_data and image parameters stored in complex_images
acq_model = AcquisitionModel(preprocessed_data, image_data);
acq_model.set_coil_sensitivity_maps(csms)

% use the acquisition model (forward projection) to simulate acquisitions
simulated_data = acq_model.forward(image_data);
sd_norm = simulated_data.norm();
% compute the difference between real and simulated acquisitions
diff = simulated_data - preprocessed_data * (sd_norm/pd_norm);
rel_residual = diff.norm()/sd_norm;
fprintf('---\n reconstruction residual norm (rel): %e\n', rel_residual)

% try to improve the reconstruction by the steepest descent step
grad = acq_model.backward(diff);
w = acq_model.forward(grad);
tau = (grad*grad)/(w*w); % locally optimal steepest descent step
refined_image_data = image_data - grad*tau;

image_array = image_data.as_array();
refined_image_array = refined_image_data.as_array();
title = 'Reconstructed image data (magnitude)';
mUtil.show_3D_array(abs(image_array), title, 'slice')
title = 'Refined image data (magnitude)';
mUtil.show_3D_array(abs(refined_image_array), title, 'slice')

