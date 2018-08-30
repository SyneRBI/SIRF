classdef NiftyResample < handle
% Class resampling nifti image using NiftyReg.

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
            name = 'SIRFRegNiftyResample';
        end
    end
    methods
        function self = NiftyResample(src)
            self.name = 'SIRFRegNiftyResample';
            self.handle_ = calllib('msirfreg', 'mSIRFReg_newObject', self.name);
            mUtilities.check_status(self.name, self.handle_)
        end
        function delete(self)
            if ~isempty(self.handle_)
                mUtilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        
        function set_reference_image(self, reference_image)
            %Set reference image.
            assert(isa(reference_image, 'mSIRFReg.ImageData'))
            mSIRFReg.setParameter(self.handle_, self.name, 'reference_image', reference_image, 'h')
        end
        function set_floating_image(self, floating_image)
            %Set floating image.
            assert(isa(floating_image, 'mSIRFReg.ImageData'))
            mSIRFReg.setParameter(self.handle_, self.name, 'floating_image', floating_image, 'h')
        end
        function set_transformation_matrix(self, filename)
            %Set transformation matrix.
            assert(ischar(filename))
            mSIRFReg.setParameter(self.handle_, self.name, 'transformation_matrix', filename)
        end
        function set_displacement_field(self, displacement_field)
            %Set displacement field.
            assert(isa(displacement_field, 'mSIRFReg.ImageDataDeformation'))
            mSIRFReg.setParameter(self.handle_, self.name, 'displacement_field', displacement_field, 'h')
        end
        function set_deformation_field(self, deformation_field)
            %Set deformation field.
            assert(isa(deformation_field, 'mSIRFReg.ImageDataDeformation'))
            mSIRFReg.setParameter(self.handle_, self.name, 'deformation_field', deformation_field, 'h')
        end
        function set_interpolation_type(self, type)
            %Set interpolation type. 0=nearest neighbour, 1=linear, 3=cubic, 4=sinc.
            assert(isinteger(type))
            mSIRFReg.setParameter(self.handle_, self.name, 'interpolation_type', type, 'i')
        end
        function set_interpolation_type_to_nearestneighbour(self)
            %Set interpolation type to nearest neighbour.
            mSIRFReg.setParameter(self.handle_, self.name, 'interpolation_type', 0, 'i')
        end
        function set_interpolation_type_to_linear(self)
            %Set interpolation type to linear.
            mSIRFReg.setParameter(self.handle_, self.name, 'interpolation_type', 1, 'i')
        end
        function set_interpolation_type_to_cubic_spline(self)
            %Set interpolation type to cubic spline.
            mSIRFReg.setParameter(self.handle_, self.name, 'interpolation_type', 3, 'i')
        end
        function set_interpolation_type_to_sinc(self)
            %Set interpolation type to sinc.
            mSIRFReg.setParameter(self.handle_, self.name, 'interpolation_type', 4, 'i')
        end
        function update(self)
            %Update.
            h = calllib('msirfreg', 'mSIRFReg_SIRFRegNiftyResample_update', self.handle_);
            mUtilities.check_status([self.name ':update'], h);
            mUtilities.delete(h)
        end
        function output = get_output(self)
            %Get output.
            output = mSIRFReg.ImageData();
            mUtilities.delete(output.handle_)
            output.handle_ = calllib('msirfreg', 'mSIRFReg_parameter', self.handle_, self.name, 'output');
            mUtilities.check_status([self.name ':get_output'], output.handle_)
        end
        function save_resampled_image(self, filename)
            %Save resampled image to file.
            assert(ischar(filename))
            h = calllib('msirfreg', 'mSIRFReg_SIRFRegNiftyResample_save_resampled_image', self.handle_, filename);
            mUtilities.check_status([self.name ':save_resampled_image'], h);
            mUtilities.delete(h)
        end
    end
end