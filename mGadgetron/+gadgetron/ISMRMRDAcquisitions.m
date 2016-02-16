classdef ISMRMRDAcquisitions < handle
    properties
        name_
        handle_
    end
    methods
        function self = ISMRMRDAcquisitions(filename)
            self.name_ = 'ISMRMRDAcquisitions';
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