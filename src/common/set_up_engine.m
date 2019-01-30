function alias = set_up_engine(engine)
% Imports a given engine.
% The engine can optionally be returned as an alias (actually a struct).
% e.g., eng=set_up_engine('STIR') enables opening STIRImageData with eng.ImageData
% Caveat: help(eng.ImageData) etc. will work, but help(eng) will not (because eng is just a struct).
% See also set_up_PET and set_up_MR and set_up_Reg

% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC.
% Copyright 2019 University College London.
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

    narginchk(1,1);
    try
        libload_sirf();
        eval(['libload_' lower(engine)])
    catch me
        fprintf(me.message)
        error('package %s failed to load\n', engine)
    end
    
    if nargout == 0 
        return;
    end
    
    %% For returning as an alias
    
    % Get the folder containing the engine
    path = fileparts(mfilename('fullpath'));
    path = fullfile(path, ['+m' engine]);
    
    % Loop over all classes and functions and set alias to handle
    files = dir([path '/*.m']);
    for i=1:size(files,1)
        [~,file,~] = fileparts(files(i).name);
        alias.(file) = eval(['@m' engine '.' file]);
    end
end