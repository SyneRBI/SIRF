classdef Gadget < handle
    properties
        handle_
        name_
    end
    methods
        function self = Gadget(fullname)
            i = strfind(fullname, '(');
            if isempty(i)
                name = strtrim(fullname);
                prop = [];
            else
                j = strfind(fullname(i + 1 : end), ')');
                name = strtrim(fullname(1 : i - 1));
                prop = fullname(i + 1 : i + j - 1);
            end
            %fprintf(name, prop);
            self.handle_ = calllib('mgadgetron', 'mGT_newObject', name);
            gadgetron.checkExecutionStatus(name, self.handle_);
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
            gadgetron.checkExecutionStatus(self.name_, handle)
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function set_properties(self, prop)
            handle = calllib('mgadgetron', 'mGT_setGadgetProperties', ...
                self.handle_, prop);
            gadgetron.checkExecutionStatus(self.name_, handle)
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
    end
end