classdef ImageWeightedMean4D < handle
% Class for performing weighted mean of deformation/displacement field images.

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
            name = 'SIRFRegImageWeightedMean4D';
        end
    end
    methods
        function self = ImageWeightedMean4D(src)
            self.name = 'SIRFRegImageWeightedMean4D';
            self.handle_ = calllib('msirfreg', 'mSIRFReg_newObject', self.name);
            mUtilities.check_status(self.name, self.handle_)
        end
        function delete(self)
            if ~isempty(self.handle_)
                mUtilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        
        function add_image(self, image, weight)
            %Add an image (filename or SIRFImageDataDeformation) and its corresponding weight.
            if isa(image, 'mSIRFReg.ImageData')
                h = calllib('msirfreg', 'mSIRFReg_SIRFRegImageWeightedMean4D_add_image', self.handle_, image.handle_, weight);
            elseif ischar(image)
                h = calllib('msirfreg', 'mSIRFReg_SIRFRegImageWeightedMean4D_add_image_filename', self.handle_, image, weight);
            else
                error("mSIRFReg.ImageWeightedMean4D.add_image: image must be SIRFImageDataDeformation or filename.")
            end
            mUtilities.check_status([self.name ':add_image'], h);
            mUtilities.delete(h)
        end
        function update(self)
            %Update.
            h = calllib('msirfreg', 'mSIRFReg_SIRFRegImageWeightedMean4D_update', self.handle_);
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
    end
end