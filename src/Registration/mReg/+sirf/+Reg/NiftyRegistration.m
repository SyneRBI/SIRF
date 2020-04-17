classdef (Abstract = true) NiftyRegistration < sirf.Reg.Registration
% Abstract class for NiftyReg registration classes.

% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2018-2020 University College London
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

    properties
        reference_image
    end
    methods(Static)
        function name = class_name()
            name = 'NiftyRegistration';
        end
    end
    methods
        function self = NiftyRegistration()
            self.name = 'NiftyRegistration';
            self.handle_ = [];
        end
        function delete(self)
            if ~isempty(self.handle_)
                sirf.Utilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function set_parameter_file(self, filename)
            %Sets the parameter filename.
            sirf.Reg.setParameter(self.handle_, 'NiftyRegistration', 'parameter_file', filename, 's')
        end
        function set_reference_mask(self, input)
            %Sets the reference mask.
            assert(isa(input, 'sirf.SIRF.ImageData'))
            sirf.Reg.setParameter(self.handle_, 'NiftyRegistration', 'reference_mask', input, 'h')
        end
        function set_floating_mask(self, input)
            %Sets the floating mask.
            assert(isa(input, 'sirf.SIRF.ImageData'))
            sirf.Reg.setParameter(self.handle_, 'NiftyRegistration', 'floating_mask', input, 'h')
        end
        function set_parameter(self, par, arg1, arg2)
            %Set string parameter. Check if any set methods match the method given by par.
            %If so, set the value given by arg. Convert to float/int etc., as necessary.
            %Up to 2 arguments, leave blank if unneeded. These are applied after parsing
            %the parameter file.
            if nargin < 3; arg1 = ''; end
            if nargin < 4; arg2 = ''; end
            h = calllib('mreg', 'mReg_NiftyRegistration_set_parameter', self.handle_, par, arg1, arg2);
        end
    end
end
