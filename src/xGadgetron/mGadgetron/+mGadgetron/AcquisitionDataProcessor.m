classdef AcquisitionDataProcessor < mGadgetron.GadgetChain
    % Class for a chain of gadgets that has AcquisitionData on input and output.
    properties
        file_
    end
    methods
        function self = AcquisitionDataProcessor(list)
%         Creates an acquisition processor specified by a list of gadgets.
%         list: Matlab cell array of gadget descriptions, each gadget 
%               description being a string of the form
%                 '[label:]gadget_name[(property1=value1[,...])]'
%               (square brackets embrace optional items, ... stands for etc.)
            self.name_ = 'AcquisitionDataProcessor';
            self.handle_ = calllib('mgadgetron', 'mGT_newObject',...
                'AcquisitionsProcessor');
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
        function acqs = process(self, input_data)
%         Returns the output from the chain for specified input.
%         input_data: AcquisitionData
            acqs = mGadgetron.AcquisitionData();
            acqs.handle_ = calllib...
                ('mgadgetron', 'mGT_processAcquisitions', ...
                self.handle_, input_data.handle_);
            mUtil.checkExecutionStatus(self.name_, acqs.handle_);
        end
    end
end