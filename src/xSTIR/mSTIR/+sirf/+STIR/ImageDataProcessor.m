classdef ImageDataProcessor < handle
% Class for image processor objects (filters).

% SyneRBI Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC.
% 
% This is software developed for the Collaborative Computational
% Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
% (http://www.ccpsynerbi.ac.uk/).
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
        name_
        handle_
        input_
        output_
    end
    methods
        function self = ImageDataProcessor()
            self.handle_ = [];
            self.input_ = [];
            self.output_ = [];
        end
        function delete(self)
            if ~isempty(self.handle_)
                %calllib('mutilities', 'mDeleteDataHandle', self.handle_)
                sirf.Utilities.delete(self.handle_)
            end
        end
        function apply(self, image)
%***SIRF*** Processes the specified image.
            h = calllib('mstir', 'mSTIR_applyImageDataProcessor',...
                self.handle_, image.handle_);
            sirf.Utilities.check_status('DataProcessor:apply', h)
            sirf.Utilities.delete(h)
            %calllib('mutilities', 'mDeleteDataHandle', h)
        end
        function set_input(self, input)
%***SIRF*** Sets the input data.
            %assert(isa(input, 'sirf.STIR.ImageData'))
            %assert(strcmp(input.class_name(), 'ImageData'))
            sirf.Utilities.assert_validity(input, 'ImageData')
            self.input_ = input;
        end
        function output = process(self, input)
%***SIRF*** Copies the input and processes the copy.
            if nargin > 1
                self.set_input(input)
            end
            if isempty(self.input_)
                error('ImageDataProcessor:input', 'input not set')
            end
            sirf.Utilities.assert_validity(self.input_, 'ImageData')
            output = self.input_.clone();
            self.apply(output)
            self.output_ = output;
        end
        function output = get_output(self)
%***SIRF*** Returns the processed data.
            output = self.output_;
        end
    end
end