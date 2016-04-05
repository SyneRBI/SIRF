classdef AcquisitionsContainer < gadgetron.DataContainer
    properties
        ordered_
    end
    methods
        function self = AcquisitionsContainer()
            self.handle_ = [];
            self.ordered_ = false;
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
                self.handle_ = [];
            end
        end
        function order(self)
            handle = calllib('mgadgetron', 'mGT_orderAcquisitions', ...
                self.handle_);
            gadgetron.checkExecutionStatus(self.name_, handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
            self.ordered_ = true;
        end
        function ordered = is_ordered(self)
            ordered = self.ordered_;
        end
    end
end