classdef MRIImageReader < handle
    properties
        handle_
        name_
    end
    methods
        function self = MRIImageReader()
            self.name_ = 'MRIImageReader';
            self.handle_ = calllib('mgadgetron', 'mGT_newObject', self.name_);
            gadgetron.checkExecutionStatus(self.name_, self.handle_);
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
            end
        end
    end
end