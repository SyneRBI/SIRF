classdef ImageDataProcessor < mGadgetron.GadgetChain
    properties
    end
    methods
        function self = ImageDataProcessor(list)
            self.name_ = 'ImageDataProcessor';
            self.handle_ = calllib...
                ('mgadgetron', 'mGT_newObject', 'ImagesProcessor');
            mUtil.checkExecutionStatus(self.name_, self.handle_);
            if nargin > 0
                for i = 1 : size(list, 2)
                    [label, name] = mUtil.label_and_name(list{i});
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
            mUtil.checkExecutionStatus(self.name_, images.handle_);
        end
    end
end