classdef QuadraticPrior < sirf.STIR.Prior
% Class for the prior that is a quadratic function of the image values.

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

    methods
        function self = QuadraticPrior()
            self.name = 'QuadraticPrior';
            self.handle_ = calllib('mstir', 'mSTIR_newObject', self.name);
            sirf.Utilities.check_status('QuadraticPrior:ctor', self.handle_)
        end
        function delete(self)
            sirf.Utilities.delete(self.handle_)
            self.handle_ = [];
        end
    end
end