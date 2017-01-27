classdef AcquisitionsProcessor < mGadgetron.GadgetChain
    properties
        file_
    end
    methods
        function self = AcquisitionsProcessor(list)
            self.name_ = 'AcquisitionsProcessor';
            self.handle_ = calllib('mgadgetron', 'mGT_newObject',...
                'AcquisitionsProcessor');
            mGadgetron.checkExecutionStatus(self.name_, self.handle_);
            if nargin > 0
                for i = 1 : size(list, 2)
                    [label, name] = mGadgetron.label_and_name(list{i});
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
%            acqs = mGadgetron.AcquisitionsContainer();
            acqs = mGadgetron.AcquisitionData();
            acqs.handle_ = calllib...
                ('mgadgetron', 'mGT_processAcquisitions', ...
                self.handle_, input_data.handle_);
            mGadgetron.checkExecutionStatus(self.name_, acqs.handle_);
        end
    end
end