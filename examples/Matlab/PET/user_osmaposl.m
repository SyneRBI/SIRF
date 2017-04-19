function user_osmaposl(engine)
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
import_str = setup_PET(engine);
eval(import_str)

try
    % direct all information printing to info.txt;
    % warning and error messages to go to Matlab Command Window
    msg_red = MessageRedirector('info.txt');

    % PET acquisition data to be read from this file
    [filename, pathname] = uigetfile('*.hs', 'Select raw data file', pet_data_path);
    acq_data = AcquisitionData(fullfile(pathname, filename));

    % create filter that zeroes the image outside a cylinder of the same
    % diameter as the image xy-section size
    filter = TruncateToCylinderProcessor();

    % create initial image estimate of dimensions and voxel sizes
    % compatible with the scanner geometry (included in the AcquisitionData
    % object ad) and initialize each voxel to 1.0
    image = acq_data.create_uniform_image(1.0);

    % apply the filter to the image
    filter.apply(image)

    % create acquisition model
    acq_model = AcquisitionModelUsingRayTracingMatrix();
    
    % create prior
    prior = QuadraticPrior();
    prior.set_penalisation_factor(0.5);

    num_subsets = 12;
    % create objective function of Poisson logarithmic likelihood type
    % compatible with the acquisition data type
    obj_fun = make_Poisson_loglikelihood(acq_data);
    obj_fun.set_acquisition_model(acq_model)
    obj_fun.set_num_subsets(num_subsets)
    obj_fun.set_up(image)
    
    num_subiterations = 2;
    image = my_osmaposl(image, obj_fun, prior, filter, num_subsets, num_subiterations);
    
    % display the reconstructed image
    image.show()
    
    image.write('my_image.hv')

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
end