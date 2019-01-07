function acquisition_sensitivity_from_attenuation(engine)
% Acquisition sensitivity model using attenuation image.

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
import_str = set_up_PET(engine);
eval(import_str)

try
    % direct all information printing to info.txt;
    % warning and error messages to go to Matlab Command Window
    MessageRedirector('info.txt', 'warn.txt');

    % raw data selected by the user is used as a template
    [filename, pathname] = uigetfile...
        ('*.hs', 'Select raw data file to be used as a template', pet_data_path);
    template = AcquisitionData(fullfile(pathname, filename));

    % create uniform acquisition data from template
    acq_data = AcquisitionData(template);
    acq_data.fill(1.0)

    % read attenuation image
    [filename, pathname] = uigetfile...
        ('*.hv', 'Select attenuation data file', pet_data_path);
    attn_image = ImageData(fullfile(pathname, filename));
    attn_image_as_array = attn_image.as_array();
    ai_dim = size(attn_image_as_array);
    z = uint16(ai_dim(3)/2);
    mUtilities.show_2D_array(attn_image_as_array(:,:,z), ...
        'Attenuation image', 'x', 'y');

    % create acquisition model
    am = AcquisitionModelUsingRayTracingMatrix();
    am.set_up(template, attn_image);

    % create acquisition sensitivity model from attenuation image
    fprintf('creating acquisition sensitivity model...\n')
    asm = AcquisitionSensitivityModel(attn_image, am);
    asm.set_up(template);
    am.set_acquisition_sensitivity(asm);

    % apply attenuation to the uniform acquisition data to obtain
    % 'bin efficiencies'
    fprintf('applying attenuation (please wait, may take a while)...\n')
    asm.unnormalise(acq_data)

    % display bin efficiencies
    acq_array = acq_data.as_array();    
    acq_dim = size(acq_array);
    z = uint16(acq_dim(3)/2);
    mUtilities.show_2D_array(acq_array(:,:,z), ...
        'Bin efficiencies', 'tang. pos.', 'views');

catch err
    % display error information
    fprintf('??? %s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
end