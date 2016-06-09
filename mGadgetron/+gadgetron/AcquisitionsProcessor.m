classdef AcquisitionsProcessor < gadgetron.GadgetChain
    properties
        file_
    end
    methods
        function self = AcquisitionsProcessor(list)
            self.name_ = 'AcquisitionsProcessor';
            self.handle_ = calllib('mgadgetron', 'mGT_newObject',...
                'AcquisitionsProcessor');
            gadgetron.checkExecutionStatus(self.name_, self.handle_);
            if nargin > 0
                for i = 1 : size(list, 2)
                    self.add_gadget(['g' num2str(i)], gadgetron.Gadget(list{i}));
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
            acqs = gadgetron.AcquisitionsContainer();
            acqs.handle_ = calllib...
                ('mgadgetron', 'mGT_processAcquisitions', ...
                self.handle_, input_data.handle_);
            gadgetron.checkExecutionStatus(self.name_, acqs.handle_);
        end
    end
end