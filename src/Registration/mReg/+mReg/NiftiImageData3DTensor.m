classdef NiftiImageData3DTensor < mReg.NiftiImageData
% Class for tensor image data.

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
            name = 'NiftiImageData3DTensor';
        end
        function obj = same_object()
            obj = mReg.NiftiImageData3DTensor();
        end
    end
    methods
        function self = NiftiImageData3DTensor(src1, src2, src3)
            narginchk(0,3)
            self.name = 'NiftiImageData3DTensor';
            if nargin < 1
                self.handle_ = calllib('mreg', 'mReg_newObject', self.name);
            elseif ischar(src1)
                self.handle_ = calllib('mreg', 'mReg_objectFromFile', self.name, src1);
            elseif nargin == 3 && isa(src1, 'mSIRF.ImageData') && isa(src2, 'mSIRF.ImageData') && isa(src3, 'mSIRF.ImageData')
                self.handle_ = calllib('mreg', 'mReg_NiftiImageData3DTensor_construct_from_3_components', self.name, src1.handle_, src2.handle_, src3.handle_);
            end
            mUtilities.check_status(self.name, self.handle_)
        end
        function delete(self)
            if ~isempty(self.handle_)
                mUtilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function write_split_xyz_components(self, filename, datatype)
            %Save to file. See nifti1.h for datatypes (e.g., float (NIFTI_TYPE_FLOAT32) = 16). 
            %Image's original datatpye is used by default.
            if nargin < 3
                datatype = -1;
            end
            h = calllib('mreg', 'mReg_NiftiImageData3DTensor_write_split_xyz_components', self.handle_, filename, datatype);
            mUtilities.check_status([self.name ':write'], h);
            mUtilities.delete(h)
        end
        function create_from_3D_image(self, src)
            %Create deformation/displacement field from 3D image.
            assert(isa(src, 'mSIRF.ImageData'), [self.name ':create_from_3D_imageInput. Input should be mSIRF.ImageData.'])
            h = calllib('mreg', 'mReg_NiftiImageData3DTensor_create_from_3D_image', self.handle_, src.handle_);
            mUtilities.check_status([self.name ':create_from_3d_image'], h);
            mUtilities.delete(h)
        end
        function flip_component(self, dim)
            % Flip component of nu
            assert(dim>=0 && dim<=2, 'Dimension to flip should be between 0 and 2.')
            h = calllib('mreg', 'mReg_NiftiImageData3DTensor_flip_component', self.handle_, dim);
            mUtilities.check_status([self.name ':flip_component'], h);
            mUtilities.delete(h)
        end
    end
end