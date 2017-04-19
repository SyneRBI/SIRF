function using_acquisition_data(engine)
% USING_ACQUISITION_DATA Demo illustrating acquisitions pre-processing 
% and displaying.
%
% In MATLAB, there are also ISMRMRD tools available for examining 
% data before processing.
%
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

% acquisitions will be read from an HDF file
[filename, pathname] = uigetfile('*.h5', 'Select raw data file', mr_data_path);
acq_data = AcquisitionData(fullfile(pathname, filename));

na = acq_data.number();
fprintf('%d acquisitions (readouts) found\n', na)

fprintf('sorting acquisitions...\n')
acq_data.sort()

% dimensions method returns size of all (i.e. including noise data) if 
% argument is passed in or if 'all' is passed in. Passing in anything else
% means not all !!
[ns, nc, na] = acq_data.dimensions('not all');

% clone acquisition data
cloned_acq_data = acq_data.clone();

% pre-process acquisition data
fprintf('processing acquisitions...\n')
processed_data = preprocess_acquisition_data(acq_data);
processed_data.sort()

% selected methods for getting information
flags = acq_data.get_info('flags');
encode_step_1 = acq_data.get_info('encode_step_1');
slice = acq_data.get_info('slice');
repetition = acq_data.get_info('repetition');

while true
    num = input('enter acquisition number (0 to stop this loop): ');
    if num < 1 || num > na
        break
    end
    fprintf('flags: %d\n', flags(num))
    fprintf('encode step 1: %d\n', encode_step_1(num))
    fprintf('slice: %d\n', slice(num))
    fprintf('repetition: %d\n', repetition(num))
end

% Data returned as complex array
acq_array0 = acq_data.as_array();

is = ns/2;
ic = nc/2;
ia = na/2;
fprintf('Value of one array element: %f\n', acq_array0(is, ic, ia))

acq_array0(is, ic, ia) = acq_array0(is, ic, ia)*10;

% Data can be replaced using fill method
acq_data.fill(acq_array0);

acq_array = acq_data.as_array();
cloned_acq_array = cloned_acq_data.as_array();
processed_array = processed_data.as_array();

fprintf('Value of same array element after replacement with 10x data: %f\n', ...
    acq_array(is, ic, ia))

acq_array = permute(acq_array, [1 3 2]);
cloned_acq_array = permute(cloned_acq_array, [1 3 2]);
processed_array = permute(processed_array, [1 3 2]);
title = 'Acquisition data (magnitude)';
mUtil.show_3D_array(abs(acq_array), title, 'samples', 'measurements', 'coil');
title = 'Cloned acquisition data (magnitude)';
mUtil.show_3D_array(abs(cloned_acq_array), title, ...
    'samples', 'measurements', 'coil');
title = 'Processed acquisition data (magnitude)';
mUtil.show_3D_array(abs(processed_array), title, ...
    'kx', 'sorted measurements', 'coil');

end


