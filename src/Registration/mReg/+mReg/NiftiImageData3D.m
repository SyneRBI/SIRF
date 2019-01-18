classdef NiftiImageData3D < mReg.NiftiImageData
% Class for 3D image data.

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
            name = 'NiftiImageData3D';
        end
    end
    methods
        function self = NiftiImageData3D(src)
            narginchk(0,1)
            self.name = 'NiftiImageData3D';
            if nargin < 1
                self.handle_ = calllib('mreg', 'mReg_newObject', self.name);
            elseif ischar(src)
                self.handle_ = calllib('mreg', 'mReg_objectFromFile', self.name, src);
            elseif isa(src, 'mSTIR.ImageData')
                self.handle_ = calllib('mreg', 'mReg_NiftiImageData3D_from_PETImageData', src.handle_);
            else
                error('NiftiImageData3D accepts no args, filename or mSTIR.ImageData.')
            end
            mUtilities.check_status(self.name, self.handle_)
        end
        function delete(self)
            if ~isempty(self.handle_)
                mUtilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function copy_data_to(self, pet_image)
            %Fill the STIRImageData with the values from NiftiImageData3D.
            assert(isa(pet_image, 'mSTIR.ImageData'))
            h = calllib('mreg', 'mReg_NiftiImageData3D_copy_data_to', self.handle_, pet_image.handle_);
            mUtilities.check_status([self.name ':copy_data_to'], h);
            mUtilities.delete(h)            
        end
    end
end