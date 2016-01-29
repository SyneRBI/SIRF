classdef GadgetChain < handle
    properties
        handle_
        name_
    end
    methods
        function self = GadgetChain()
            self.name_ = 'GadgetChain';
            self.handle_ = calllib('mgadgetron', 'mGT_newObject', self.name_);
            gadgetron.checkExecutionStatus(self.name_, self.handle_);
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
            end
        end
        function add_reader(self, id, reader)
            handle = calllib...
                ('mgadgetron', 'mGT_addReader', self.handle_, id, reader.handle_);
            gadgetron.checkExecutionStatus(self.name_, handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function add_writer(self, id, writer)
            handle = calllib...
                ('mgadgetron', 'mGT_addWriter', self.handle_, id, writer.handle_);
            gadgetron.checkExecutionStatus(self.name_, handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function add_gadget(self, id, gadget)
            handle = calllib...
                ('mgadgetron', 'mGT_addGadget', self.handle_, id, gadget.handle_);
            gadgetron.checkExecutionStatus(self.name_, handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
    end
end