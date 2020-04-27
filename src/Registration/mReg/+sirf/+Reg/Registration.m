classdef (Abstract = true) Registration < handle
% Abstract class for registration classes.

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
        name
        handle_
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
                sirf.Utilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function set_reference_image(self, input)
            %Sets the reference image.
            assert(isa(input, 'sirf.SIRF.ImageData'))
            self.reference_image = input;
            sirf.Reg.setParameter(self.handle_, 'Registration', 'reference_image', input, 'h')
        end
        function set_floating_image(self, input)
            %Sets the floating image. Will clear any previous floating images.
            assert(isa(input, 'sirf.SIRF.ImageData'))
            sirf.Reg.setParameter(self.handle_, 'Registration', 'floating_image', input, 'h')
        end
        function add_floating_image(self, input)
            %Add floating image.
            assert(isa(input, 'sirf.SIRF.ImageData'))
            h = calllib('mreg', 'mReg_Registration_add_floating', self.handle_, input.handle_);
            sirf.Utilities.check_status([self.name ':add_floating_image'], h);
            sirf.Utilities.delete(h)
        end
        function set_reference_image_filename(self, filename)
            %Sets the reference image filename.
            assert(ischar(filename))
            self.reference_image = sirf.Reg.NiftiImageData(filename);
            h = calllib('mreg', 'mReg_Registration_set_reference_image_filename', self.handle_, filename);
            sirf.Utilities.check_status([self.name ':set_reference_image_filename'], h);
            sirf.Utilities.delete(h)
        end
        function set_floating_image_filename(self, filename)
            %Sets the floating image filename.
            assert(ischar(filename))
            h = calllib('mreg', 'mReg_Registration_set_floating_image_filename', self.handle_, filename);
            sirf.Utilities.check_status([self.name ':set_floating_image_filename'], h);
            sirf.Utilities.delete(h)
        end
        function add_floating_image_filename(self, filename)
            %Add floating image filename.
            assert(ischar(filename))
            h = calllib('mreg', 'mReg_Registration_add_floating_image_filename', self.handle_, filename);
            sirf.Utilities.check_status([self.name ':add_floating_image_filename'], h);
            sirf.Utilities.delete(h)
        end
        function clear_floating_images(self)
            %Clear floating images.
            h = calllib('mreg', 'mReg_Registration_clear_floatings', self.handle_);
            sirf.Utilities.check_status([self.name ':clear_floating_images'], h);
            sirf.Utilities.delete(h)
        end
        function output = get_output(self, idx)
            %Gets the registered image. 1-based.
            assert(~isempty(self.reference_image) && ~isempty(self.reference_image.handle_))
            if nargin < 2; idx=1; end
            output = self.reference_image.same_object();
            sirf.Utilities.delete(output.handle_)
            output.handle_ = calllib('mreg', 'mReg_Registration_get_output', self.handle_, round(idx-1));
            sirf.Utilities.check_status([self.name ':get_output'], output.handle_)
        end
        function process(self)
            %Run the registration.
            assert(~isempty(self.handle_), 'Registration.process: Registration object is empty.')
            h = calllib('mreg', 'mReg_Registration_process', self.handle_);
            sirf.Utilities.check_status([self.name ':process'], h);
            sirf.Utilities.delete(h)
        end
        function output = get_deformation_field_forward(self, idx)
            %Gets the forward deformation field image. 1-based.
            if nargin < 2; idx=1; end
            output = sirf.Reg.NiftiImageData3DDeformation();
            output.handle_ = calllib('mreg', 'mReg_Registration_get_deformation_displacement_image', self.handle_, 'forward_deformation', round(idx-1));
            sirf.Utilities.check_status([self.name ':get_deformation_field_forward'], output.handle_);
        end
        function output = get_deformation_field_inverse(self, idx)
            %Gets the inverse deformation field image. 1-based.
            if nargin < 2; idx=1; end
            output = sirf.Reg.NiftiImageData3DDeformation();
            output.handle_ = calllib('mreg', 'mReg_Registration_get_deformation_displacement_image', self.handle_, 'inverse_deformation', round(idx-1));
            sirf.Utilities.check_status([self.name ':get_deformation_field_inverse'], output.handle_);
        end
        function output = get_displacement_field_forward(self, idx)
            %Gets the forward displacement field image. 1-based.
            if nargin < 2; idx=1; end
            output = sirf.Reg.NiftiImageData3DDisplacement();
            output.handle_ = calllib('mreg', 'mReg_Registration_get_deformation_displacement_image', self.handle_, 'forward_displacement', round(idx-1));
            sirf.Utilities.check_status([self.name ':get_displacement_field_forward'], output.handle_);
        end
        function output = get_displacement_field_inverse(self, idx)
            %Gets the inverse displacement field image. 1-based.
            if nargin < 2; idx=1; end
            output = sirf.Reg.NiftiImageData3DDisplacement();
            output.handle_ = calllib('mreg', 'mReg_Registration_get_deformation_displacement_image', self.handle_, 'inverse_displacement', round(idx-1));
            sirf.Utilities.check_status([self.name ':get_displacement_field_inverse'], output.handle_);
        end
    end
end
