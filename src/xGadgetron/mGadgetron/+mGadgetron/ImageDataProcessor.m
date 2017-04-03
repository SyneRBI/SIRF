classdef ImageDataProcessor < mGadgetron.GadgetChain
% Class for a chain of gadgets that has ImageData on input and output.

% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC.
% Copyright 2015 - 2017 University College London.
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
    end
    methods
        function self = ImageDataProcessor(list)
%         Creates an image processor specified by a list of gadgets.
%         list: Matlab cell array of gadget descriptions, each gadget 
%               description being a string of the form
%                 '[label:]gadget_name[(property1=value1[,...])]'
%               (square brackets embrace optional items, ... stands for etc.)
            self.name_ = 'ImageDataProcessor';
            self.handle_ = calllib...
                ('mgadgetron', 'mGT_newObject', 'ImagesProcessor');
            mUtil.checkExecutionStatus(self.name_, self.handle_);
            if nargin > 0
                for i = 1 : size(list, 2)
                    [label, name] = mUtil.label_and_name(list{i});
                    self.add_gadget(label, mGadgetron.Gadget(name));
                end
            end
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
            end
            self.handle_ = [];
        end
        function images = process(self, input_data)
%         Returns the output from the chain for specified input.
%         input_data: ImageData
            images = mGadgetron.ImageData();
            images.handle_ = calllib...
                ('mgadgetron', 'mGT_processImages', ...
                self.handle_, input_data.handle_);
            mUtil.checkExecutionStatus(self.name_, images.handle_);
        end
    end
end