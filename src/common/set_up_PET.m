function alias = set_up_PET(engine)
% Imports a given PET engine. If no engine is given, a default is used.
% The engine can optionally be returned as an alias (actually a struct).
% e.g., eng=set_up_engine('STIR') enables opening STIRImageData with eng.ImageData
% Caveat: help(eng.ImageData) etc. will work, but help(eng) will not (because eng is just a struct).
% See also set_up_engine and set_up_MR and set_up_Reg

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

    if nargin == 0 || isempty(engine)
        engine = 'STIR';
    end
    if nargout == 0
    	set_up_engine(engine);
    else
    	alias = set_up_engine(engine);
    end
end