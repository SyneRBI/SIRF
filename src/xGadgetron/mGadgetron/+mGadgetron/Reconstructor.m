classdef Reconstructor < mGadgetron.GadgetChain
% Class for reconstructing images using Gadgetron 
% 
% Reconstructor Methods:
%    Reconstructor - can accept a list of gadgets for gadgetron chain
%    set_input(input_data)  - sets the input data for recon
%    process  - performs the call to gadgetron
%    get_output  - returns output

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
        input_
        images_
    end
    methods
        function self = Reconstructor(list)
%         Creates an image processor chain defined by an optional argument,
%         a Matlab cell array of gadget descriptions, each description 
%         being a Matlab string of the form
%             '[label:]gadget_name[(property1=value1[,...])]'
%         (square brackets embrace optional items, ... stands for etc.)
%         If no argument is present, the empty chain is created.
%         The use of labels enables subsequent setting of gadget properties 
%         using set_gadget_property(label, property, value).
            self.name_ = 'ImagesReconstructor';
            self.handle_ = calllib('mgadgetron', 'mGT_newObject', self.name_);
            self.input_ = [];
            self.images_ = [];
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
        function set_input(self, input_data)
%***SIRF*** Sets the specified AcquisitionData argument as the input.
%         See also PROCESS
            self.input_ = input_data;
        end
        function process(self)
%***SIRF*** Calls Gadgetron.
            if isempty(self.input_)
                error('MRIReconstruction:no_input', ...
                    'no input data for reconstruction')
            end
            self.images_ = mGadgetron.ImageData();
            self.images_.handle_ = calllib...
                ('mgadgetron', 'mGT_reconstructImages', ...
                self.handle_, self.input_.handle_);
            mUtil.checkExecutionStatus(self.name_, self.images_.handle_);
        end
        function images = get_output(self, subset)
%***SIRF*** get_output(subset) returns the results of the image reconstruction 
%         as an ImageData object;
%         subset: either 'image' or 'gfactor'
            images = self.images_;
            if nargin > 1
                images = images.select('GADGETRON_DataRole', subset);
            end
        end
    end
end
