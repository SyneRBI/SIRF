classdef (Abstract = true) Transformation < handle & matlab.mixin.Heterogeneous
% Abstract class for transformations.

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

    properties
        name
        handle_
    end
    methods(Static)
        function name = class_name()
            name = 'SIRFRegTransformation';
        end
    end
    methods
        function self = Transformation()
            self.name = 'SIRFRegTransformation';
            self.handle_ = [];
        end
        function delete(self)
            if ~isempty(self.handle_)
                mUtilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function output = get_as_deformation_field(self, ref)
            %Get any type of transformation as a deformation field.
            %This is useful for joining them together. Require a reference
            %image for converting transformation matrices to deformations.
            assert(isa(ref, 'mSIRFReg.NiftiImage3D'))
            output = mSIRFReg.NiftiImage3DDeformation();
            output.handle_ = calllib('msirfreg', 'mSIRFReg_SIRFRegTransformation_get_as_deformation_field', self.handle_, ref.handle_);
            mUtilities.check_status([self.name ':get_as_deformation_field'], output.handle_);
        end
    end
end