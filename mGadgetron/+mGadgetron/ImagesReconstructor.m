classdef ImagesReconstructor < mGadgetron.GadgetChain
% ImagesReconstructor 
% Class for reconstructing images using Gadgetron 
% 
% ImagesReconstructor Methods:
%    ImagesReconstructor - can accept a list of gadgets for gadgetron chain
%    set_input(input_data)  - sets the input data for recon
%    process  - performs the call to gadgetron
%    get_output  - returns output

    properties
        input_
        images_
    end
    methods
        function self = ImagesReconstructor(list)
        % Accepts a cell array of gadget names as input to form the Gadgetron chain.
        % Each element of list is a string of the form 
        % 'LABEL:gadget_name' with the 'LABEL:' optional. Use of a label
        % enables subsequent setting of gadget property using set_gadget_property
            self.name_ = 'ImagesReconstructor';
            self.handle_ = calllib('mgadgetron', 'mGT_newObject', self.name_);
            self.input_ = [];
            self.images_ = [];
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
        function set_input(self, input_data)
         % set_input - Sets the Acquisition Data used for the recon
         % See also PROCESS
            self.input_ = input_data;
        end
        function process(self)
        % process - Calls the Gadgetron
            if isempty(self.input_)
                error('MRIReconstruction:no_input', ...
                    'no input data for reconstruction')
            end
            self.images_ = mGadgetron.ImagesContainer();
            self.images_.handle_ = calllib...
                ('mgadgetron', 'mGT_reconstructImages', ...
                self.handle_, self.input_.handle_);
            mGadgetron.checkExecutionStatus(self.name_, self.images_.handle_);
        end
        function images = get_output(self, subset)
        % get_output - Returns an ImagesContainer?
            images = self.images_;
            if nargin > 1
                images = images.select('GADGETRON_DataRole', subset);
            end
        end
    end
end