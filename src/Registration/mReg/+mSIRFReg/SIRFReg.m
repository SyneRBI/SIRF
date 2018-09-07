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
            assert(isa(input, 'mSIRFReg.ImageData'))
            mSIRFReg.setParameter(self.handle_, 'SIRFReg', 'reference_image', input, 'h')
        end
        function set_floating_image(self, input)
            %Sets the floating image.
            assert(isa(input, 'mSIRFReg.ImageData'))
            mSIRFReg.setParameter(self.handle_, 'SIRFReg', 'floating_image', input, 'h')
        end
        function output = get_output(self)
            %Gets the registered image.
            output = mSIRFReg.ImageData();
            mUtilities.delete(output.handle_)
            output.handle_ = calllib('msirfreg', 'mSIRFReg_parameter', self.handle_, 'SIRFReg', 'output');
            mUtilities.check_status([self.name ':get_output'], output.handle_)
        end
        function save_warped_image(self, filename)
            %Save warped image.
            h = calllib('msirfreg', 'mSIRFReg_SIRFReg_save_image', self.handle_, filename);
            mUtilities.check_status([self.name ':save_warped_image'], h);
            mUtilities.delete(h)
        end
        function save_deformation_field_fwrd(self, filename, split_xyz)
            %Save forward deformation field image to file.
            h = calllib('msirfreg', 'mSIRFReg_SIRFReg_save_deformation_displacement_image', self.handle_, filename, 'fwrd_deformation', split_xyz);
            mUtilities.check_status([self.name ':save_deformation_field_fwrd'], h);
            mUtilities.delete(h)
        end
        function save_deformation_field_back(self, filename, split_xyz)
            %Save backward deformation field image to file.
            h = calllib('msirfreg', 'mSIRFReg_SIRFReg_save_deformation_displacement_image', self.handle_, filename, 'back_deformation', split_xyz);
            mUtilities.check_status([self.name ':save_deformation_field_back'], h);
            mUtilities.delete(h)
        end
        function save_displacement_field_fwrd(self, filename, split_xyz)
            %Save forward displacement field image to file.
            h = calllib('msirfreg', 'mSIRFReg_SIRFReg_save_deformation_displacement_image', self.handle_, filename, 'fwrd_displacement', split_xyz);
            mUtilities.check_status([self.name ':save_displacement_field_fwrd'], h);
            mUtilities.delete(h)
        end
        function save_displacement_field_back(self, filename, split_xyz)
            %Save backward displacement field image to file.
            h = calllib('msirfreg', 'mSIRFReg_SIRFReg_save_deformation_displacement_image', self.handle_, filename, 'back_displacement', split_xyz);
            mUtilities.check_status([self.name ':save_displacement_field_back'], h);
            mUtilities.delete(h)
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
            output = mSIRFReg.ImageDataDeformation();
            output.handle_ = calllib('msirfreg', 'mSIRFReg_SIRFReg_get_deformation_displacement_image', self.handle_, 'fwrd_deformation');
            mUtilities.check_status([self.name ':get_deformation_field_fwrd'], output.handle_);
        end
        function output = get_deformation_field_back(self)
            %Gets the backwards deformation field image.
            output = mSIRFReg.ImageDataDeformation();
            output.handle_ = calllib('msirfreg', 'mSIRFReg_SIRFReg_get_deformation_displacement_image', self.handle_, 'back_deformation');
            mUtilities.check_status([self.name ':get_deformation_field_back'], output.handle_);
        end
        function output = get_displacement_field_fwrd(self)
            %Gets the forward displacement field image.
            output = mSIRFReg.ImageDataDeformation();
            output.handle_ = calllib('msirfreg', 'mSIRFReg_SIRFReg_get_deformation_displacement_image', self.handle_, 'fwrd_displacement');
            mUtilities.check_status([self.name ':get_displacement_field_fwrd'], output.handle_);
        end
        function output = get_displacement_field_back(self)
            %Gets the backwards displacement field image.
            output = mSIRFReg.ImageDataDeformation();
            output.handle_ = calllib('msirfreg', 'mSIRFReg_SIRFReg_get_deformation_displacement_image', self.handle_, 'back_displacement');
            mUtilities.check_status([self.name ':get_displacement_field_back'], output.handle_);
        end
    end
end