classdef AcquisitionsProcessor < gadgetron.GadgetChain
    properties
        file_
    end
    methods
        function self = AcquisitionsProcessor()
            self.name_ = 'AcquisitionsProcessor';
            self.file_ = [datestr(now,'dd_mm_yyyy_HH_MM_SS') '.h5'];
            self.handle_ = calllib('mgadgetron', ...
                'mGT_acquisitionsProcessor', self.file_);
            gadgetron.checkExecutionStatus(self.name_, self.handle_);
        end
        function delete(self)
            calllib('mutilities', 'mDeleteObject', self.handle_)
            self.handle_ = [];
        end
        function acqs = process(self, input_data)
            acqs = gadgetron.ISMRMRDAcquisitions();
            acqs.handle_ = calllib...
                ('mgadgetron', 'mGT_processAcquisitions', ...
                self.handle_, input_data.handle_);
            gadgetron.checkExecutionStatus(self.name_, acqs.handle_);
        end
    end
end