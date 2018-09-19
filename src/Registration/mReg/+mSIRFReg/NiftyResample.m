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
            assert(isa(reference_image, 'mSIRFReg.NiftiImage3D'), 'NiftyResample::set_reference_image expects NiftiImage3D')
            mSIRFReg.setParameter(self.handle_, self.name, 'reference_image', reference_image, 'h')
        end
        function set_floating_image(self, floating_image)
            %Set floating image.
            assert(isa(floating_image, 'mSIRFReg.NiftiImage3D'), 'NiftyResample::set_floating_image expects NiftiImage3D')
            mSIRFReg.setParameter(self.handle_, self.name, 'floating_image', floating_image, 'h')
        end
        function add_transformation_affine(self, src)
            %Set transformation matrix.
            assert(isa(src, 'mSIRFReg.TransformationAffine'), 'NiftyResample::add_transformation_affine expects TransformationAffine.')
            h = calllib('msirfreg', 'mSIRFReg_SIRFRegNiftyResample_add_transformation', self.handle_, src.handle_, 'affine');
        end

        function add_transformation_disp(self, src)
            %Set displacement field.
            assert(isa(src, 'mSIRFReg.TransformationDisplacement'), 'NiftyResample::add_transformation_disp expects TransformationDisplacement.')
            h = calllib('msirfreg', 'mSIRFReg_SIRFRegNiftyResample_add_transformation', self.handle_, src.handle_, 'displacement');
        end

        function add_transformation_def(self, src)
            %Set deformation field.
            assert(isa(src, 'mSIRFReg.TransformationDeformation'), 'NiftyResample::add_transformation_def expects TransformationDeformation.')
            h = calllib('msirfreg', 'mSIRFReg_SIRFRegNiftyResample_add_transformation', self.handle_, src.handle_, 'deformation');
        end
        function set_interpolation_type(self, type)
            %Set interpolation type. 0=nearest neighbour, 1=linear, 3=cubic, 4=sinc.
            mSIRFReg.setParameter(self.handle_, self.name, 'interpolation_type', type, 'i')
        end
        function set_interpolation_type_to_nearest_neighbour(self)
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
            output = mSIRFReg.NiftiImage3D();
            mUtilities.delete(output.handle_)
            output.handle_ = calllib('msirfreg', 'mSIRFReg_parameter', self.handle_, self.name, 'output');
            mUtilities.check_status([self.name ':get_output'], output.handle_)
        end
    end
end
