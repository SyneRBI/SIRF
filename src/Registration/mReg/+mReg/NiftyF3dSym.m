classdef NiftyF3dSym < mReg.Registration
% Registration class using NiftyReg's symmetric f3d.

% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2018-2019 University College London
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
            mUtilities.check_status(self.name, self.handle_)
        end
        function delete(self)
            if ~isempty(self.handle_)
                mUtilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function set_floating_time_point(self, floating_time_point)
            %Set floating time point.
            mReg.setParameter(self.handle_, self.name, 'floating_time_point', floating_time_point, 'i')
        end
        function set_reference_time_point(self, reference_time_point)
            %Set reference time point.
            mReg.setParameter(self.handle_, self.name, 'reference_time_point', reference_time_point, 'i')
        end
        function set_initial_affine_transformation(self, src)
            %Set initial affine transformation.
            assert(isa(src, 'mReg.AffineTransformation'))
            mReg.setParameter(self.handle_, self.name, 'initial_affine_transformation', src, 'h');
        end
    end
end