function acquisition_data_from_scanner_info(engine)
% A demo showing basics of PET acquisition data handling.

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

% all acquisition data generated by this script will be stored in memory
% (the input data remain in the input file);
% default storage scheme 'file' keeps all acquisition data generated by
% the script in scratch files deleted after the script terminates
AD = PET.AcquisitionData();
scheme = AD.get_storage_scheme();
AD.set_storage_scheme('memory');
% scheme = AcquisitionData.get_storage_scheme();
% AcquisitionData.set_storage_scheme('memory');

try
    % create acquisition data from scanner parameters
    acq_data = PET.AcquisitionData('Siemens_mMR');
    % set all values to 1.0
    acq_data.fill(1.0);

    % copy the acquisition data into a Python array
    acq_array = acq_data.as_array();
    acq_dim = size(acq_array);
    fprintf('acquisition data dimensions (maximum resolution): %d x %d x %d\n', acq_dim)

    % create acquisition data from scanner parameters but with axial compression etc
    span=11;
    max_ring_diff=-1;
    view_mash_factor=2;
    acq_data = PET.AcquisitionData('Siemens_mMR', span, max_ring_diff, view_mash_factor);
    % copy the acquisition data into a Python array
    acq_array = acq_data.as_array();
    acq_dim = size(acq_array);
    fprintf('acquisition data dimensions (span 11, view mashing 2): %d x %d x %d\n', acq_dim)

    % write the acquisition data to a file (commented out for this demo)
    % acq_data.write('example_mMR_ones.hs')

catch err
    % display error information
    fprintf('??? %s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
AD.set_storage_scheme(scheme);
%AcquisitionData.set_storage_scheme(scheme);