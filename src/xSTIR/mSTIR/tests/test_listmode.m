function [failed, ntests] = test_listmode(record, engine)
% PET test set 1.

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
if nargin < 2
    engine = [];
end
if nargin < 1
    record = false;
end

pet = set_up_PET(engine);

% define raw data source
data_path = sirf.Utilities.examples_data_path('PET');
raw_data_file = fullfile(data_path, 'mMR', 'list.l.hdr');

lm2sino = pet.ListmodeToSinograms();
lm2sino.set_input(raw_data_file);

prompt_rate_threshold = 73036.;
known_time = 22.;

time_at_which_prompt_rate_exceeds_threshold = ...
    lm2sino.get_time_at_which_prompt_rate_exceeds_threshold(prompt_rate_threshold);

assert(abs(time_at_which_prompt_rate_exceeds_threshold-known_time) <= 1e-4, ...
	'ListmodeToSinograms::get_time_at_which_prompt_rate_exceeds_threshold failed')

failed = 0;
ntests = 1;

end
