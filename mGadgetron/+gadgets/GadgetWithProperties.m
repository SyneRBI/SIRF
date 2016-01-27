classdef GadgetWithProperties < handle
    properties
        name_
    end
    methods
        function set_property(self, prop, value)
            handle = calllib('mgadgetron', 'mGT_setGadgetProperty', ...
                self.handle_, prop, value);
            gadgetron.checkExecutionStatus(self.name_, handle)
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
    end
end