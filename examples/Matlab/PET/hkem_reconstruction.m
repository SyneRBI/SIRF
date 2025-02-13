function hkem_reconstruction(engine)
% HKEM reconstruction demo. 

% SyneRBI Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2021 Rutherford Appleton Laboratory STFC.
% Copyright 2015 - 2019 University College London.
% 
% This is software developed for the Collaborative Computational
% Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
% (http://www.ccpsynerbi.ac.uk/).
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

if nargin < 1
    engine = [];
end
% import_str = set_up_PET(engine);
% eval(import_str)
PET = set_up_PET(engine);
AD = PET.AcquisitionData();
AD.set_storage_scheme('memory');
%AcquisitionData.set_storage_scheme('memory');
pet_data_path = sirf.Utilities.examples_data_path('PET');

try
    % direct all printing to MatlabCommand Window
    PET.MessageRedirector('stdout');

    % create acquisition model
    acq_model = PET.AcquisitionModelUsingRayTracingMatrix();
    
    % PET acquisition data to be read from this file
    [filename, pathname] = uigetfile('*.hs', 'Select raw data file', pet_data_path);
    acq_data = PET.AcquisitionData(fullfile(pathname, filename));

    % read anatomical image
    [filename, pathname] = uigetfile('*.hv', 'Select anatomical image file', ...
        pet_data_path);
    anatomical_image = PET.ImageData(fullfile(pathname, filename));
    ai_arr = anatomical_image.as_array();
    ai_arr(ai_arr < 0) = 0;
    anatomical_image.fill(ai_arr);
    fprintf('anatomical image dimensions:%d %d %d\n', ...
        anatomical_image.dimensions());
    anatomical_image.show();

    % create initial image estimate
    image = anatomical_image.get_uniform_copy();

    fprintf('setting up acquisition model...\n')
    acq_model.set_up(acq_data, image);

    % create objective function of Poisson logarithmic likelihood type
    % compatible with the acquisition data type
    obj_fun = PET.make_Poisson_loglikelihood(acq_data);
    obj_fun.set_acquisition_model(acq_model)
    
    num_subsets = 4;
    num_subiterations = 2;
    
    % select Kernelized Ordered Subsets Maximum A-Posteriori One Step Late
    % as the reconstruction algorithm
    recon = PET.KOSMAPOSLReconstructor();
    recon.set_objective_function(obj_fun);
    recon.set_num_subsets(num_subsets);
    recon.set_input(acq_data);
    recon.set_anatomical_prior(anatomical_image);
    recon.set_num_neighbours(5);
    recon.set_num_non_zero_features(3);
    recon.set_sigma_m(2.0);
    recon.set_sigma_p(3.0);
    recon.set_sigma_dm(4.0);
    recon.set_sigma_dp(5.0);
    recon.set_only_2D(true);
    recon.set_hybrid(true);

    fprintf('setting up, please wait...\n')
    recon.set_up(image)

    for iter = 1 : num_subiterations
        fprintf('\n------------- Subiteration %d\n', recon.get_subiteration_num())
        % perform a sub-iteration
        recon.update(image)
        % display the current image at z = 10
        image.show(10);
    end

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
end
