function run_all_scripts(path)
% This script contains common part of all run_all.m scripts
% to reduce duplication.

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

filelist = dir(fullfile(path, '*.m'));
nf = length(filelist);
for k = 1 : nf
    filename = filelist(k).name;
    if ~strcmpi(filename, 'run_all.m')
        fprintf('Running %s\n', filename)
        try
            run(fullfile(path, filename));
            fprintf('ok\n')
        catch
            fprintf('failed\n')
        end
    end
end