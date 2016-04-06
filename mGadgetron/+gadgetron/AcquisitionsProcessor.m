classdef AcquisitionsProcessor < gadgetron.GadgetChain
    properties
        file_
    end
    methods
        function self = AcquisitionsProcessor()
            self.name_ = 'AcquisitionsProcessor';
            self.handle_ = calllib('mgadgetron', 'mGT_newObject',...
                'AcquisitionsProcessor');
            gadgetron.checkExecutionStatus(self.name_, self.handle_);
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