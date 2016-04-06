classdef AcquisitionsContainer < gadgetron.DataContainer
    properties
        sorted_
    end
    methods
        function self = AcquisitionsContainer()
            self.handle_ = [];
            self.sorted_ = false;
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
                self.handle_ = [];
            end
        end
        function sort(self)
            handle = calllib('mgadgetron', 'mGT_orderAcquisitions', ...
                self.handle_);
            gadgetron.checkExecutionStatus(self.name_, handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
            self.sorted_ = true;
        end
        function sorted = is_sorted(self)
            sorted = self.sorted_;
        end
    end
end