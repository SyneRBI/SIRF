classdef ImageFinishGadget < handle
    properties
        handle_
        name_
    end
    methods
        function self = ImageFinishGadget()
            self.name_ = 'ImageFinishGadget';
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