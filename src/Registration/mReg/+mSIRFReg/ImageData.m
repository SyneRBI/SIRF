classdef ImageData < handle
% Class for image data.

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
            name = 'SIRFImageData';
        end
    end
    methods
        function self = ImageData(src)
            narginchk(0,1)
            self.name = 'SIRFImageData';
            if nargin < 1
                self.handle_ = calllib('msirfreg', 'mSIRFReg_newObject', self.name);
            elseif ischar(src)
                self.handle_ = calllib('msirfreg', 'mSIRFReg_objectFromFile', self.name, src);
            elseif isa(src, 'mSTIR.ImageData')
                self.handle_ = calllib('msirfreg', 'mSIRFReg_SIRFImageData_from_PETImageData', src.handle_);
            else
                error('ImageData accepts no args, filename or mSTIR.ImageData.')
            end
            mUtilities.check_status(self.name, self.handle_)
        end
        function delete(self)
            if ~isempty(self.handle_)
                mUtilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function save_to_file(self, filename)
            %Save to file.
            h = calllib('msirfreg', 'mSIRFReg_SIRFImageData_save_to_file', self.handle_, filename);
            mUtilities.check_status([self.name ':save_to_file'], h);
            mUtilities.delete(h)
        end
        function value = get_max(self)
            %Get max.
            value = mSIRFReg.parameter(self.handle_, self.name, 'max', 'f');
        end
        function value = get_min(self)
            %Get min.
            value = mSIRFReg.parameter(self.handle_, self.name, 'max', 'f');
        end
        function copy_data_to(self, pet_image)
            %Fill the STIRImageData with the values from SIRFImageData.
            assert(isa(pet_image, 'mSTIR.ImageData'))
            h = calllib('msirfreg', 'mSIRFReg_SIRFImageData_copy_data_to', self.handle_, pet_image.handle_);
            mUtilities.check_status([self.name ':copy_data_to'], h);
            mUtilities.delete(h)            
        end
    end
end