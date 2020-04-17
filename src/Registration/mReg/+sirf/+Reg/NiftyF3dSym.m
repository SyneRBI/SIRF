classdef NiftyF3dSym < sirf.Reg.NiftyRegistration
% Registration class using NiftyReg's symmetric f3d.

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

    methods(Static)
        function name = class_name()
            name = 'NiftyF3dSym';
        end
    end
    methods
        function self = NiftyF3dSym()
            self.name = 'NiftyF3dSym';
            self.handle_ = calllib('mreg', 'mReg_newObject', self.name);
            sirf.Utilities.check_status(self.name, self.handle_)
        end
        function delete(self)
            if ~isempty(self.handle_)
                sirf.Utilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function set_floating_time_point(self, floating_time_point)
            %Set floating time point.
            sirf.Reg.setParameter(self.handle_, self.name, 'floating_time_point', floating_time_point, 'i')
        end
        function set_reference_time_point(self, reference_time_point)
            %Set reference time point.
            sirf.Reg.setParameter(self.handle_, self.name, 'reference_time_point', reference_time_point, 'i')
        end
        function set_initial_affine_transformation(self, src)
            %Set initial affine transformation.
            assert(isa(src, 'sirf.Reg.AffineTransformation'))
            sirf.Reg.setParameter(self.handle_, self.name, 'initial_affine_transformation', src, 'h');
        end
    end
    methods(Static)
	function print_all_wrapped_methods()
	    %Print list of all wrapped methods.
	    disp('In C++, this class is templated. "dataType" corresponds to "float" for Matlab and python.')
	    h = calllib('mreg', 'mReg_NiftyRegistration_print_all_wrapped_methods', 'NiftyF3dSym');
	    sirf.Utilities.check_status('parameter', h)
	end
    end
end