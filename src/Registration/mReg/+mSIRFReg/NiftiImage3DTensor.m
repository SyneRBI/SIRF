classdef NiftiImage3DTensor < mSIRFReg.NiftiImage
% Class for tensor image data.

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

    methods(Static)
        function name = class_name()
            name = 'NiftiImage3DTensor';
        end
    end
    methods
        function self = NiftiImage3DTensor(src1, src2, src3)
            narginchk(0,3)
            self.name = 'NiftiImage3DTensor';
            if nargin < 1
                self.handle_ = calllib('msirfreg', 'mSIRFReg_newObject', self.name);
            elseif ischar(src1)
                self.handle_ = calllib('msirfreg', 'mSIRFReg_objectFromFile', self.name, src1);
            elseif nargin == 3 && isa(src1, 'mSIRFReg.NiftiImage3D') && isa(src2, 'mSIRFReg.NiftiImage3D') && isa(src3, 'mSIRFReg.NiftiImage3D')
                self.handle_ = calllib('msirfreg', 'mSIRFReg_NiftiImage3DTensor_construct_from_3_components', self.name, src1.handle_, src2.handle_, src3.handle_);
            end
            mUtilities.check_status(self.name, self.handle_)
        end
        function delete(self)
            if ~isempty(self.handle_)
                mUtilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function save_to_file_split_xyz_components(self, filename)
            % Save to file.
            h = calllib('msirfreg', 'mSIRFReg_NiftiImage3DTensor_save_to_file_split_xyz_components', self.handle_, filename);
            mUtilities.check_status([self.name ':save_to_file'], h);
            mUtilities.delete(h)
        end
        function create_from_3D_image(self, src)
            %Create deformation/displacement field from 3D image.
            assert(isa(src, 'mSIRFReg.NiftiImage3D'), [self.name ':create_from_3D_imageInput. Input should be NiftiImage3D.'])
            h = calllib('msirfreg', 'mSIRFReg_NiftiImage3DTensor_create_from_3D_image', self.handle_, src.handle_);
            mUtilities.check_status([self.name ':create_from_3d_image'], h);
            mUtilities.delete(h)
        end
    end
end