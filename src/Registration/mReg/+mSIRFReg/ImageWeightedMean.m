classdef ImageWeightedMean < handle
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
            name = 'SIRFRegImageWeightedMean';
        end
    end
    methods
        function self = ImageWeightedMean(src)
            self.name = 'SIRFRegImageWeightedMean';
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
            %Add an image (filename or SIRFImageData) and its corresponding weight.
            if isa(image, 'mSIRFReg.ImageData')
                h = calllib('msirfreg', 'mSIRFReg_SIRFRegImageWeightedMean_add_image', self.handle_, image.handle_, weight);
            elseif ischar(image)
                h = calllib('msirfreg', 'mSIRFReg_SIRFRegImageWeightedMean_add_image_filename', self.handle_, image, weight);
            else
                error("mSIRFReg.ImageWeightedMean.add_image: image must be SIRFImageData or filename.")
            end
            mUtilities.check_status([self.name ':add_image'], h);
            mUtilities.delete(h)
        end
        function update(self)
            %Update.
            h = calllib('msirfreg', 'mSIRFReg_SIRFRegImageWeightedMean_update', self.handle_);
            mUtilities.check_status([self.name ':update'], h);
            mUtilities.delete(h)
        end
        function save_image_to_file(self, filename)
            %Save image to file.
            assert(ischar(filename))
            h = calllib('msirfreg', 'mSIRFReg_SIRFRegImageWeightedMean_save_image_to_file', self.handle_, filename);
            mUtilities.check_status([self.name ':save_image_to_file'], h);
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