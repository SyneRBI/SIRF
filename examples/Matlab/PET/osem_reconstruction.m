function osem_reconstruction(engine)
% OSEM reconstruction demo. 
% We actually use the OSMAPOSL reconstructor in this demo. This reconstructor
% implements an Ordered Subsets (OS) version of the One Step Late algorithm (OSL)
% from Green et al for Maximum a Posteriori (MAP) maximisation. Here we use it
% for Maximum Likelihood (ML) in which case it is equivalent to OSEM.

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

if nargin < 1
    engine = [];
end
import_str = set_up_PET(engine);
eval(import_str)
pet_data_path = mUtilities.examples_data_path('PET');

try
    % direct all printing to MatlabCommand Window
    MessageRedirector('stdout');

    % create acquisition model
    acq_model = AcquisitionModelUsingRayTracingMatrix();
    
    % PET acquisition data to be read from this file
    [filename, pathname] = uigetfile('*.hs', 'Select raw data file', pet_data_path);
    acq_data = AcquisitionData(fullfile(pathname, filename));

    % create initial image estimate of dimensions and voxel sizes
    % compatible with the scanner geometry (included in the AcquisitionData
    % object ad) and initialize each voxel to 1.0
    image = acq_data.create_uniform_image(1.0);

    % create objective function of Poisson logarithmic likelihood type
    % compatible with the acquisition data type
    obj_fun = make_Poisson_loglikelihood(acq_data);
    obj_fun.set_acquisition_model(acq_model)
    
    num_subiterations = 2;
    
    % select Ordered Subsets Maximum A-Posteriori One Step Late as the
    % reconstruction algorithm (since we are not using a penalty, or prior, in
    % this example, we actually run OSEM);
    % this algorithm does not converge to the maximum of the objective function
    % but is used in practice to speed-up calculations
    recon = OSMAPOSLReconstructor();    
    recon.set_objective_function(obj_fun)
    recon.set_input(acq_data)
    recon.set_num_subsets(12)
    
    % set up the reconstructor based on a sample image
    % (checks the validity of parameters, sets up objective function
    % and other objects involved in the reconstruction, which involves
    % computing/reading sensitivity image etc etc.)
    fprintf('setting up, please wait...\n')
    recon.set_up(image)

    % display the initial image
    z = 20;
    image_array = image.as_array();
    mUtilities.show_2D_array(image_array(:,:,z), 'initial image', 'x', 'y');

    % set the initial image estimate
    recon.set_current_estimate(image)

    % suppress further information printing
    MessageRedirector();

    % in order to see the reconstructed image evolution
    % open up the user's access to the iterative process
    % rather than allow recon.reconstruct to do all job at once
    for iter = 1 : num_subiterations
        fprintf('\n--------------------- Subiteration %d\n', iter)
        % perform an iteration
        recon.update_current_estimate()
        % display the current image
        image_array = recon.get_current_estimate().as_array();
        the_title = sprintf('iteration %d', iter);
        mUtilities.show_2D_array(image_array(:,:,z), the_title, 'x', 'y');
    end

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
end