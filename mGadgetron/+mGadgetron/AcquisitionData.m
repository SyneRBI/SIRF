classdef AcquisitionData < mGadgetron.AcquisitionsContainer
    properties
        name_
    end
    methods
        function self = AcquisitionData(filename)
            self.name_ = 'MR_Acquisitions';
            self.handle_ = [];
            if nargin > 0
                self.handle_ = calllib('mgadgetron', ...
                    'mGT_ISMRMRDAcquisitionsFromFile', filename);
                mGadgetron.checkExecutionStatus(self.name_, self.handle_);
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