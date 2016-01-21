classdef SimpleReconGadget < handle
    properties
        handle_
        name_
    end
    methods
        function self = SimpleReconGadget()
            self.name_ = 'SimpleReconGadget';
            self.handle_ = calllib('mgadgetron', 'mNewObject', self.name_);
            gadgetron.checkExecutionStatus(self.name_, self.handle_);
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mgadgetron', 'mDeleteObject', self.handle_)
            end
        end
    end
end