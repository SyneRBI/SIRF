classdef ImagesReconstructor < gadgetron.GadgetChain
    properties
        input_
        images_
    end
    methods
        function self = ImagesReconstructor()
            self.name_ = 'ImageReconstructor';
            self.handle_ = calllib('mgadgetron', 'mGT_newObject', self.name_);
            self.input_ = [];
            self.images_ = [];
            gadgetron.checkExecutionStatus(self.name_, self.handle_);
        end
        function delete(self)
            calllib('mutilities', 'mDeleteObject', self.handle_)
            self.handle_ = [];
        end
        function set_input(self, input_data)
            self.input_ = input_data;
        end
        function process(self)
            if isempty(self.input_)
                error('MRIReconstruction:no_input', ...
                    'no input data for reconstruction')
            end
            self.images_ = gadgetron.ImagesContainer();
            self.images_.handle_ = calllib...
                ('mgadgetron', 'mGT_reconstructImages', ...
                self.handle_, self.input_.handle_);
            gadgetron.checkExecutionStatus(self.name_, self.images_.handle_);
        end
        function images = get_output(self)
            images = self.images_;
        end
    end
end