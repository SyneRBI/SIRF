function import_str = set_up_Reg(engine, alias)
% Creates a string, evaluating which imports PET engine.
% Optionally creates also its alias (e.g. named PET, to have PET.ImageData etc.)

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

    if isempty(engine)
        engine = 'SIRFReg';
    end
    try
        eval(['libload_' lower(engine)])
    catch me
        fprintf(me.message)
        error('package %s failed to load\n', engine)
    end
    if nargin < 2
        import_str = ['import m' engine '.*'];
    else
        if ~strcmp(['m' engine], alias)
            filename = mfilename();
            filepath = mfilename('fullpath');
            l = length(filepath) - length(filename);
            path = filepath(1:l);
            copyfile([path '/+m' engine], [path '/+' alias], 'f')
        end
        import_str = ' ';
    end
end