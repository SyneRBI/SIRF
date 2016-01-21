classdef ISMRMRDataset < handle
    properties
        name_
        handle_
        header_
    end
    methods
        function self = ISMRMRDataset(file)
            self.name_ = 'ISMRMRDataset';
            self.handle_ = calllib...
                ('mgadgetron', 'mGT_ISMRMRDatasetFromFile', file, '/dataset');
            gadgetron.checkExecutionStatus(self.name_, self.handle_);
            self.header_ = calllib('mgadgetron', 'mNewObject', 'string');
            handle = calllib...
                ('mgadgetron', 'mGT_readISMRMRDatasetHeader', ...
                self.handle_, self.header_);
            gadgetron.checkExecutionStatus(self.name_, handle);
            calllib('mgadgetron', 'mDeleteDataHandle', handle)
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mgadgetron', 'mDeleteObject', self.handle_)
            end
        end
        function read_header(self)
            handle = calllib('mgadgetron', 'mGT_readISMRMRDatasetHeader', ...
                self.handle_, self.header_);
            gadgetron.checkExecutionStatus(self.name_, handle);
            calllib('mgadgetron', 'mDeleteDataHandle', handle)
        end
    end
end