classdef AcquisitionDataProcessor < sirf.Gadgetron.GadgetChain
% Class for a chain of gadgets that has AcquisitionData on input and output.

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
        input_
        output_
    end
    methods
        function self = AcquisitionDataProcessor(list)
%         Creates an acquisition processor chain defined by an optional argument,
%         a Matlab cell array of gadget descriptions, each description 
%         being a Matlab string of the form
%             '[label:]gadget_name[(property1=value1[,...])]'
%         (square brackets embrace optional items, ... stands for etc.)
%         If no argument is present, the empty chain is created.
%         The use of labels enables subsequent setting of gadget properties 
%         using set_gadget_property(label, property, value).
            self.name_ = 'AcquisitionDataProcessor';
            self.handle_ = calllib('mgadgetron', 'mGT_newObject',...
                'AcquisitionsProcessor');
            sirf.Utilities.check_status(self.name_, self.handle_);
            if nargin > 0
                for i = 1 : size(list, 2)
                    [label, name] = sirf.Utilities.label_and_name(list{i});
                    self.add_gadget(label, sirf.Gadgetron.Gadget(name));
                end
            end
        end
        function delete(self)
            if ~isempty(self.handle_)
                sirf.Utilities.delete(self.handle_)
                %calllib('mutilities', 'mDeleteObject', self.handle_)
            end
            self.handle_ = [];
        end
        function set_input(self, input)
%***SIRF*** Sets the input data.
            %assert(isa(input, 'sirf.Gadgetron.AcquisitionData'))
            %assert(strcmp(input.class_name(), 'AcquisitionData'))
            sirf.Utilities.assert_validity(input, 'AcquisitionData')
            self.input_ = input;
        end
        function acqs = process(self, input_data)
%***SIRF*** Returns the output from the chain 
%         for the input specified by the agrument.
%         both input and output are of type AcquisitionData.
            if nargin > 1 % input is supplied as an argument
                self.set_input(input_data)
            end
            if isempty(self.input_)
                error('ImageDataProcessor:input', 'input not set')
            end
            sirf.Utilities.assert_validity(self.input_, 'AcquisitionData')
            acqs = sirf.Gadgetron.AcquisitionData();
            acqs.handle_ = calllib...
                ('mgadgetron', 'mGT_processAcquisitions', ...
                self.handle_, self.input_.handle_);
            sirf.Utilities.check_status(self.name_, acqs.handle_);
            self.output_ = acqs;
        end
        function output = get_output(self)
%***SIRF*** Returns the processed data.
            output = self.output_;
        end
    end
end