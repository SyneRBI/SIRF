classdef MR_Acquisitions < gadgetron.AcquisitionsContainer
    properties
        name_
    end
    methods
        function self = MR_Acquisitions(filename)
            self.name_ = 'MR_Acquisitions';
            self.handle_ = [];
            if nargin > 0
                self.handle_ = calllib('mgadgetron', ...
                    'mGT_ISMRMRDAcquisitionsFromFile', filename);
                gadgetron.checkExecutionStatus(self.name_, self.handle_);
            end
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
                self.handle_ = [];
            end
        end
    end
end