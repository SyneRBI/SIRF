function using_acquisition_data(engine)
% A demo showing basics of PET acquisition data handling.

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
    % select acquisition data to test
    [filename, pathname] = uigetfile...
        ('*.hs', 'Select raw data file', pet_data_path);
    acq_data = AcquisitionData(fullfile(pathname, filename));

    % copy the acquisition data into a Matlab array
    acq_array = acq_data.as_array();

    scale = max(max(max(acq_array)))/4;
    acq_dim = size(acq_array);
    x = acq_dim(1)/2;
    y = acq_dim(2)/2;
    z = acq_dim(3)/2;

    % display the acquisition data
    figure
    imshow(acq_array(:,:,z)/scale);
    title('acquisition data')

    % clone the acquisition data
    new_acq_data = acq_data.clone();
    % display the cloned data
    acq_array = new_acq_data.as_array();
    figure
    imshow(acq_array(:,:,z)/scale);
    title('acquisition data cloned')

    % fill the cloned data with the acquisition data multiplied by 10
    % and see the difference at (x, y, z)
    fprintf('acq_data at (%d,%d,%d): %f\n', x, y, z, acq_array(x, y, z))
    new_acq_data.fill(10*acq_array)
    acq_array = new_acq_data.as_array();
    fprintf('new_acq_data at (%d,%d,%d): %f\n', x, y, z, acq_array(x, y, z))
    
catch err
    % display error information
    fprintf('??? %s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
end