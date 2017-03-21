classdef Gadget < handle
    % Class for Gadgetron gadgets.
    properties
        handle_
        name_
    end
    methods
        function self = Gadget(fullname)
%         Creates a gadget of specified type and properties.
%         name: a string of the form gadget_type(property1=value1, ...)
            [name, prop] = mUtil.name_and_parameters(fullname);
            self.handle_ = calllib('mgadgetron', 'mGT_newObject', name);
            mUtil.checkExecutionStatus(name, self.handle_);
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
%         Assigns specified value to specified gadget property.
%         prop : property name (string)
%         value: property value (string)
            handle = calllib('mgadgetron', 'mGT_setGadgetProperty', ...
                self.handle_, prop, value);
            mUtil.checkExecutionStatus(self.name_, handle)
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function set_properties(self, prop)
%         Assigns specified values to specified gadget properties.
%         prop: a string with comma-separated list of property value assignments 
%               prop_name=prop_value
            handle = calllib('mgadgetron', 'mGT_setGadgetProperties', ...
                self.handle_, prop);
            mUtil.checkExecutionStatus(self.name_, handle)
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
    end
end