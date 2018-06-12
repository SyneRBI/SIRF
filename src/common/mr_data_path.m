function p = mr_data_path
% Tries to find path to MR raw data.
% The user may like to set a Matlab variable SIRF_MR_DATA_PATH
% to the path to their raw MR data.
% If it is not set, the path to SIRF subfolder /data/examples/MR
% will be used.

% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC.
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

SIRF_MR_DATA_PATH = getenv('SIRF_MR_DATA_PATH');
if ~isempty(SIRF_MR_DATA_PATH)
    p = SIRF_MR_DATA_PATH;
else
    SIRF_PATH = getenv('SIRF_PATH');
    if ~isempty(SIRF_PATH)
        p = [SIRF_PATH '/data/examples/MR'];
    else
        p = './';
    end
end
end