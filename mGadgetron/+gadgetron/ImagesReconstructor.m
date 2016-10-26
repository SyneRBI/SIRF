classdef ImagesReconstructor < gadgetron.GadgetChain
    properties
        input_
        images_
    end
    methods
        function self = ImagesReconstructor(list)
            self.name_ = 'ImagesReconstructor';
            self.handle_ = calllib('mgadgetron', 'mGT_newObject', self.name_);
            self.input_ = [];
            self.images_ = [];
            gadgetron.checkExecutionStatus(self.name_, self.handle_);
            if nargin > 0
                for i = 1 : size(list, 2)
                    [label, name] = gadgetron.label_and_name(list{i});
                    self.add_gadget(label, gadgetron.Gadget(name));
                end
            end
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
            end
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