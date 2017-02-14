classdef Gadget < handle
    properties
        handle_
        name_
    end
    methods
        function self = Gadget(fullname)
            [name, prop] = mGadgetron.name_and_parameters(fullname);
            self.handle_ = calllib('mgadgetron', 'mGT_newObject', name);
            mGadgetron.checkExecutionStatus(name, self.handle_);
            if ~isempty(prop)
                self.set_properties(prop)
            end
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
            mGadgetron.checkExecutionStatus(self.name_, handle)
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function set_properties(self, prop)
            handle = calllib('mgadgetron', 'mGT_setGadgetProperties', ...
                self.handle_, prop);
            mGadgetron.checkExecutionStatus(self.name_, handle)
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
    end
end