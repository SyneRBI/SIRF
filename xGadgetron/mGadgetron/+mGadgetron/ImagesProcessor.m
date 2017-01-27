classdef ImagesProcessor < mGadgetron.GadgetChain
    properties
    end
    methods
        function self = ImagesProcessor(list)
            self.name_ = 'ImagesProcessor';
            self.handle_ = calllib('mgadgetron', 'mGT_newObject', self.name_);
            mGadgetron.checkExecutionStatus(self.name_, self.handle_);
            if nargin > 0
                for i = 1 : size(list, 2)
                    [label, name] = mGadgetron.label_and_name(list{i});
                    self.add_gadget(label, mGadgetron.Gadget(name));
                end
            end
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
            end
            self.handle_ = [];
        end
        function images = process(self, input_data)
            images = mGadgetron.ImageData();
            images.handle_ = calllib...
                ('mgadgetron', 'mGT_processImages', ...
                self.handle_, input_data.handle_);
            mGadgetron.checkExecutionStatus(self.name_, images.handle_);
        end
    end
end