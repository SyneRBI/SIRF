classdef NiftiImageData3DDisplacement < mReg.NiftiImageData3DTensor & mReg.Transformation
% Class for displacement image data.

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
            name = 'NiftiImageData3DDisplacement';
        end
    end
    methods
        function self = NiftiImageData3DDisplacement(src1, src2, src3)
            narginchk(0,3)
            self.name = 'NiftiImageData3DDisplacement';
            if nargin < 1
                self.handle_ = calllib('msirfreg', 'mReg_newObject', self.name);
            elseif ischar(src1)
                self.handle_ = calllib('msirfreg', 'mReg_objectFromFile', self.name, src1);
            elseif nargin == 3 && isa(src1, 'mReg.NiftiImageData3D') && isa(src2, 'mReg.NiftiImageData3D') && isa(src3, 'mReg.NiftiImageData3D')
                self.handle_ = calllib('msirfreg', 'mReg_NiftiImageData3DTensor_construct_from_3_components', self.name, src1.handle_, src2.handle_, src3.handle_);
            elseif isa(src1, 'mReg.NiftiImageData3DDeformation')
                self.handle_ = calllib('msirfreg', 'mReg_NiftiImageData3DDisplacement_create_from_def', src1.handle_);
            end
            mUtilities.check_status(self.name, self.handle_)
        end
        function delete(self)
            if ~isempty(self.handle_)
                mUtilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
    end
end