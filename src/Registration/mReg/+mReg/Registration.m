classdef (Abstract = true) Registration < handle
% Abstract class for registration classes.

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

    properties
        name
        handle_
        reference_image
    end
    methods(Static)
        function name = class_name()
            name = 'Registration';
        end
    end
    methods
        function self = Registration()
            self.name = 'Registration';
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
            mReg.setParameter(self.handle_, 'Registration', 'parameter_file', filename, 's')
        end
        function set_reference_image(self, input)
            %Sets the reference image.
            assert(isa(input, 'mSIRF.ImageData'))
            self.reference_image = input;
            mReg.setParameter(self.handle_, 'Registration', 'reference_image', input, 'h')
        end
        function set_floating_image(self, input)
            %Sets the floating image.
            assert(isa(input, 'mSIRF.ImageData'))
            mReg.setParameter(self.handle_, 'Registration', 'floating_image', input, 'h')
        end
        function set_reference_mask(self, input)
            %Sets the reference mask.
            assert(isa(input, 'mSIRF.ImageData'))
            mReg.setParameter(self.handle_, 'Registration', 'reference_mask', input, 'h')
        end
        function set_floating_mask(self, input)
            %Sets the floating mask.
            assert(isa(input, 'mSIRF.ImageData'))
            mReg.setParameter(self.handle_, 'Registration', 'floating_mask', input, 'h')
        end
        function output = get_output(self)
            %Gets the registered image.
            assert(~isempty(self.reference_image) && ~isempty(self.reference_image.handle_))
            output = self.reference_image.same_object();
            mUtilities.delete(output.handle_)
            output.handle_ = calllib('mreg', 'mReg_parameter', self.handle_, 'Registration', 'output');
            mUtilities.check_status([self.name ':get_output'], output.handle_)
        end
        function process(self)
            %Run the registration.
            assert(~isempty(self.handle_), 'Registration.process: Registration object is empty.')
            h = calllib('mreg', 'mReg_Registration_process', self.handle_);
            mUtilities.check_status([self.name ':process'], h);
            mUtilities.delete(h)
        end
        function output = get_deformation_field_forward(self)
            %Gets the forward deformation field image.
            output = mReg.NiftiImageData3DDeformation();
            output.handle_ = calllib('mreg', 'mReg_Registration_get_deformation_displacement_image', self.handle_, 'forward_deformation');
            mUtilities.check_status([self.name ':get_deformation_field_forward'], output.handle_);
        end
        function output = get_deformation_field_inverse(self)
            %Gets the inverse deformation field image.
            output = mReg.NiftiImageData3DDeformation();
            output.handle_ = calllib('mreg', 'mReg_Registration_get_deformation_displacement_image', self.handle_, 'inverse_deformation');
            mUtilities.check_status([self.name ':get_deformation_field_inverse'], output.handle_);
        end
        function output = get_displacement_field_forward(self)
            %Gets the forward displacement field image.
            output = mReg.NiftiImageData3DDisplacement();
            output.handle_ = calllib('mreg', 'mReg_Registration_get_deformation_displacement_image', self.handle_, 'forward_displacement');
            mUtilities.check_status([self.name ':get_displacement_field_forward'], output.handle_);
        end
        function output = get_displacement_field_inverse(self)
            %Gets the inverse displacement field image.
            output = mReg.NiftiImageData3DDisplacement();
            output.handle_ = calllib('mreg', 'mReg_Registration_get_deformation_displacement_image', self.handle_, 'inverse_displacement');
            mUtilities.check_status([self.name ':get_displacement_field_inverse'], output.handle_);
        end
        function set_parameter(self, par, arg1, arg2)
            %Set string parameter. Check if any set methods match the method given by par.
            %If so, set the value given by arg. Convert to float/int etc., as necessary.
            %Up to 2 arguments, leave blank if unneeded. These are applied after parsing
            %the parameter file.
            if nargin < 3; arg1 = ''; end
            if nargin < 4; arg2 = ''; end
            h = calllib('mreg', 'mReg_Registration_set_parameter', self.handle_, par, arg1, arg2);
        end
    end
end
