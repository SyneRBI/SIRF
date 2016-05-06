classdef Gadget < handle
    properties
        handle_
        name_
    end
    methods
        function self = Gadget(name)
            self.handle_ = calllib('mgadgetron', 'mGT_newObject', name);
            gadgetron.checkExecutionStatus(name, self.handle_);
            self.name_ = name;
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
            end
        end
        function set_property(self, prop, value)
            handle = calllib('mgadgetron', 'mGT_setGadgetProperty', ...
                self.handle_, prop, value);
            gadgetron.checkExecutionStatus(self.name_, handle)
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
    end
end