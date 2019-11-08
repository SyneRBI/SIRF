function listmode_to_sinograms(engine)
% Listmode-to-sinograms conversion demo.
% Notes:
%   output will be in the current working directory.
%   input defaults to the mMR subfolder of pet_data_path

% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2018 Rutherford Appleton Laboratory STFC.
% Copyright 2018 University College London.
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
pet_data_path = sirf.Utilities.examples_data_path('PET');

AD = PET.AcquisitionData();
AD.set_storage_scheme('memory');
%AcquisitionData.set_storage_scheme('memory');

try
    % direct all information printing to info.txt;
    % warning and error messages to go to Matlab Command Window
    PET.MessageRedirector('info.txt', 'warn.txt');

    % create listmode-to-sinograms converter object
    lm2sino = PET.ListmodeToSinograms();

    default_path = fullfile(pet_data_path, 'mMR');

    [filename, pathname] = uigetfile...
        ('*.l.hdr', 'Select listmode data file', default_path);
    list_file = fullfile(pathname, filename);
    % get the filename of a template AcquisitionData
    % the template is used to specify the sizes of the output sinogram.
    % see the acquisition_data_from_scanner_info demo for an example how to 
    % make your own template file
    [filename, pathname] = uigetfile...
        ('*.hs', 'Select raw data file to be used as a template', default_path);
    tmpl_file = fullfile(pathname, filename);
    
    % set input, output and template files
    lm2sino.set_input(list_file)
    lm2sino.set_output_prefix('sinograms')
    lm2sino.set_template(tmpl_file)

    % set interval
    lm2sino.set_time_interval(0, 10)

    % set some flags as examples (the following values are the defaults)
    lm2sino.flag_on('store_prompts')
    lm2sino.flag_off('interactive')

    % set up the converter
    lm2sino.set_up()

    % convert
    lm2sino.process()

    % get access to the sinograms
    acq_data = lm2sino.get_output();
    % copy the acquisition data into a Python array
    acq_array = acq_data.as_array();
    %acq_dim = size(acq_array);
    acq_dim = acq_data.dimensions();
    fprintf('acquisition data dimensions: %d x %d x %d x %d\n', acq_dim)
    z = uint16(acq_dim(3)/2);
    sirf.Utilities.show_2D_array(acq_array(:,:,z), ...
        'acquisition data', 'tang. pos.', 'views');

    % compute randoms
    fprintf('estimating randoms, please wait...\n')
    randoms = lm2sino.estimate_randoms();
    rnd_array = randoms.as_array();
    sirf.Utilities.show_2D_array(rnd_array(:,:,z), ...
        'randoms', 'tang. pos.', 'views');

catch err
    % display error information
    fprintf('??? %s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
