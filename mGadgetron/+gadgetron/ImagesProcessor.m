classdef ImagesProcessor < gadgetron.GadgetChain
    properties
    end
    methods
        function self = ImagesProcessor()
            self.name_ = 'ImagesProcessor';
            self.handle_ = calllib('mgadgetron', 'mGT_newObject', self.name_);
            gadgetron.checkExecutionStatus(self.name_, self.handle_);
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
            end
            self.handle_ = [];
        end
        function images = process(self, input_data)
            images = gadgetron.ImagesContainer();
            images.handle_ = calllib...
                ('mgadgetron', 'mGT_processImages', ...
                self.handle_, input_data.handle_);
            gadgetron.checkExecutionStatus(self.name_, images.handle_);
        end
    end
end