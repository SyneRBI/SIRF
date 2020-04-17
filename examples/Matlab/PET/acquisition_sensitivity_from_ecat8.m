function acquisition_sensitivity_from_ecat8(engine)
% Acquisition sensitivity model using ECAT8 bin normalization.

% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2018 Rutherford Appleton Laboratory STFC.
% Copyright 2015 - 2018 University College London.
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
% import_str = set_up_PET(engine);
% eval(import_str)
PET = set_up_PET(engine);
pet_data_path = [sirf.Utilities.examples_data_path('PET') '/mMR'];

try
    % direct all information printing to info.txt;
    % warning and error messages to go to Matlab Command Window
    PET.MessageRedirector('info.txt', 'warn.txt');

    % raw data selected by the user is used as a template
    [filename, pathname] = uigetfile...
        ('*.hs', 'Select raw data file to be used as a template', pet_data_path);
    template = PET.AcquisitionData(fullfile(pathname, filename));

    % create acquisition sensitivity model from ECAT8 normalization data
    [filename, pathname] = uigetfile...
        ('*.n.hdr', 'Select ECAT8 normalization file', pet_data_path);
    asm = PET.AcquisitionSensitivityModel(fullfile(pathname, filename));
    asm.set_up(template);

    % create a uniform acquisition data from template
    acq_data = PET.AcquisitionData(template);
    acq_data.fill(1.0)

    % apply normalization to the uniform acquisition data to obtain
    % bin efficiencies
    fwd_data = asm.forward(acq_data);

    % display bin efficiencies
    acq_array = fwd_data.as_array();    
    acq_dim = size(acq_array);
    z = uint16(acq_dim(3)/2);
    sirf.Utilities.show_2D_array(acq_array(:,:,z), ...
        'Bin efficiencies', 'tang. pos.', 'views');

catch err
    % display error information
    fprintf('??? %s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
end