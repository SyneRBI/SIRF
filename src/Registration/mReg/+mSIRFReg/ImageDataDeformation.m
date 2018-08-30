classdef ImageDataDeformation < mSIRFReg.ImageData
% Class for deformation/displacement image data.

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
            name = 'SIRFImageDataDeformation';
        end
    end
    methods
        function self = ImageDataDeformation(filename)
            narginchk(0,1)
            self.name = 'SIRFImageDataDeformation';
            if nargin < 1
                self.handle_ = calllib('msirfreg', 'mSIRFReg_newObject', self.name);
            else
                self.handle_ = calllib('msirfreg', 'mSIRFReg_objectFromFile', self.name, filename);
            end
            mUtilities.check_status(self.name, self.handle_)
        end
        function delete(self)
            if ~isempty(self.handle_)
                mUtilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function save_to_file(self, filename, split_xyz)
            % Save to file.
            h = calllib('msirfreg', 'mSIRFReg_SIRFImageDataDeformation_save_to_File', self.handle_, filename);
            mUtilities.check_status([self.name ':save_to_file'], h);
            mUtilities.delete(h)
        end
        function create_from_3D_image(self, src)
            %Create deformation/displacement field from 3D image.
            assert(isa(src, 'mSIRFReg.ImageData'))
            h = calllib('msirfreg', 'mSIRFReg_SIRFImageDataDeformation_create_from_3D_image', self.handle_, src.handle_);
            mUtilities.check_status([self.name ':create_from_3D_image'], h);
            mUtilities.delete(h)
        end
    end
end