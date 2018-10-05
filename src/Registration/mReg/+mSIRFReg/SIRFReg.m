classdef (Abstract = true) SIRFReg < handle
% Abstract class for registration classes.

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
            name = 'SIRFReg';
        end
    end
    methods
        function self = SIRFReg()
            self.name = 'SIRFReg';
            self.handle_ = [];
        end
        function delete(self)
            if ~isempty(self.handle_)
                mUtilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function set_parameter_file(self, filename)
            %Sets the parameter filename.
            mSIRFReg.setParameter(self.handle_, 'SIRFReg', 'parameter_file', filename, 's')
        end
        function set_reference_image(self, input)
            %Sets the reference image.
            assert(isa(input, 'mSIRFReg.NiftiImage3D'))
            mSIRFReg.setParameter(self.handle_, 'SIRFReg', 'reference_image', input, 'h')
        end
        function set_floating_image(self, input)
            %Sets the floating image.
            assert(isa(input, 'mSIRFReg.NiftiImage3D'))
            mSIRFReg.setParameter(self.handle_, 'SIRFReg', 'floating_image', input, 'h')
        end
        function output = get_output(self)
            %Gets the registered image.
            output = mSIRFReg.NiftiImage3D();
            mUtilities.delete(output.handle_)
            output.handle_ = calllib('msirfreg', 'mSIRFReg_parameter', self.handle_, 'SIRFReg', 'output');
            mUtilities.check_status([self.name ':get_output'], output.handle_)
        end
        function update(self)
            %Run the registration.
            assert(~isempty(self.handle_), 'SIRFReg.update: Registration object is empty.')
            h = calllib('msirfreg', 'mSIRFReg_SIRFReg_update', self.handle_);
            mUtilities.check_status([self.name ':update'], h);
            mUtilities.delete(h)
        end
        function output = get_deformation_field_fwrd(self)
            %Gets the forward deformation field image.
            output = mSIRFReg.NiftiImage3DDeformation();
            output.handle_ = calllib('msirfreg', 'mSIRFReg_SIRFReg_get_deformation_displacement_image', self.handle_, 'fwrd_deformation');
            mUtilities.check_status([self.name ':get_deformation_field_fwrd'], output.handle_);
        end
        function output = get_deformation_field_back(self)
            %Gets the backwards deformation field image.
            output = mSIRFReg.NiftiImage3DDeformation();
            output.handle_ = calllib('msirfreg', 'mSIRFReg_SIRFReg_get_deformation_displacement_image', self.handle_, 'back_deformation');
            mUtilities.check_status([self.name ':get_deformation_field_back'], output.handle_);
        end
        function output = get_displacement_field_fwrd(self)
            %Gets the forward displacement field image.
            output = mSIRFReg.NiftiImage3DDisplacement();
            output.handle_ = calllib('msirfreg', 'mSIRFReg_SIRFReg_get_deformation_displacement_image', self.handle_, 'fwrd_displacement');
            mUtilities.check_status([self.name ':get_displacement_field_fwrd'], output.handle_);
        end
        function output = get_displacement_field_back(self)
            %Gets the backwards displacement field image.
            output = mSIRFReg.NiftiImage3DDisplacement();
            output.handle_ = calllib('msirfreg', 'mSIRFReg_SIRFReg_get_deformation_displacement_image', self.handle_, 'back_displacement');
            mUtilities.check_status([self.name ':get_displacement_field_back'], output.handle_);
        end
        function set_parameter(self, par, arg1, arg2)
            %Set string parameter. Check if any set methods match the method given by par.
            %If so, set the value given by arg. Convert to float/int etc., as necessary.
            %Up to 2 arguments, leave blank if unneeded. These are applied after parsing
            %the parameter file.
            if nargin < 3; arg1 = ''; end
            if nargin < 4; arg2 = ''; end
            h = calllib('msirfreg', 'mSIRFReg_SIRFReg_set_parameter', self.handle_, par, arg1, arg2);
        end
    end
end