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
        function set_gadget_property(self, id, prop, value)
            if ischar(value)
                v = value;
            elseif islogical(value)
                if value
                    v = 'true';
                else
                    v = 'false';
                end
            else
                v = num2str(value);
            end
            hg = calllib('mgadgetron', 'mGT_parameter', ...
                self.handle_, 'gadget_chain', id);
            gadgetron.checkExecutionStatus(self.name_, hg);
            hv = calllib('mgadgetron', 'mGT_setGadgetProperty', ...
                hg, prop, v);
            gadgetron.checkExecutionStatus(self.name_, hv)
            calllib('mutilities', 'mDeleteDataHandle', hg)
            calllib('mutilities', 'mDeleteDataHandle', hv)
        end
        function v = value_of_gadget_property(self, id, prop)
            hg = calllib('mgadgetron', 'mGT_parameter', ...
                self.handle_, 'gadget_chain', id);
            gadgetron.checkExecutionStatus(self.name_, hg);
            hv = calllib('mgadgetron', 'mGT_parameter', hg, 'gadget', prop);
            gadgetron.checkExecutionStatus(self.name_, hv);
            v = calllib('mutilities', 'mCharDataFromHandle', hv);
            calllib('mutilities', 'mDeleteDataHandle', hg)
            calllib('mutilities', 'mDeleteDataHandle', hv)
        end
    end
end