classdef MRIReconstruction < gadgetron.GadgetChain
    properties
    end
    methods
        function self = MRIReconstruction()
            self.name_ = 'MRIReconstruction';
            self.handle_ = calllib('mgadgetron', 'mGT_newObject', self.name_);
            gadgetron.checkExecutionStatus(self.name_, self.handle_);
        end
        function delete(self)
            calllib('mutilities', 'mDeleteObject', self.handle_)
            self.handle_ = [];
        end
        function process(self, input_data)
            handle = calllib...
                ('mgadgetron', 'mGT_runMRIReconstruction', ...
                self.handle_, input_data.handle_);
            gadgetron.checkExecutionStatus(self.name_, handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function images = get_output(self)
            imgs.handle_ = calllib...
                ('mgadgetron', 'mGT_reconstructedImagesList', self.handle_);
            gadgetron.checkExecutionStatus(self.name_, imgs.handle_);
            images = gadgetron.ImagesList(imgs);
        end
    end
end