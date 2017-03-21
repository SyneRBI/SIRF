classdef GadgetChain < handle
    % Class for Gadgetron chains.
    properties
        handle_
        name_
    end
    methods
        function self = GadgetChain()
            self.name_ = 'GadgetChain';
            self.handle_ = calllib('mgadgetron', 'mGT_newObject', self.name_);
            mUtil.checkExecutionStatus(self.name_, self.handle_);
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
            end
        end
        function add_reader(self, id, reader)
%         Adds reader gadget (a gadget that receives data from the client) to the
%         chain.
%         id    : gadget id (string)
%         reader: Gadget of reader type
            handle = calllib...
                ('mgadgetron', 'mGT_addReader', self.handle_, id, reader.handle_);
            mUtil.checkExecutionStatus(self.name_, handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function add_writer(self, id, writer)
%         Adds writer gadget (a gadget that sends data to the client) to the
%         chain.
%         id    : gadget id (string)
%         writer: Gadget of writer type
            handle = calllib...
                ('mgadgetron', 'mGT_addWriter', self.handle_, id, writer.handle_);
            mUtil.checkExecutionStatus(self.name_, handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function add_gadget(self, id, gadget)
%         Adds a gadget to the chain.
%         id    : gadget id (string)
%         writer: Gadget
            handle = calllib...
                ('mgadgetron', 'mGT_addGadget', self.handle_, id, gadget.handle_);
            mUtil.checkExecutionStatus(self.name_, handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function set_gadget_property(self, id, prop, value)
%         Assigns specified value to specified gadget property.
%         id   : gadget id
%         prop : property name (string)
%         value: property value (string)
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
            mUtil.checkExecutionStatus(self.name_, hg);
            hv = calllib('mgadgetron', 'mGT_setGadgetProperty', ...
                hg, prop, v);
            mUtil.checkExecutionStatus(self.name_, hv)
            calllib('mutilities', 'mDeleteDataHandle', hg)
            calllib('mutilities', 'mDeleteDataHandle', hv)
        end
        function v = value_of_gadget_property(self, id, prop)
%         Returns the string representation of the value of specified property.
%         id  : gadget id
%         prop: property name (string)
            hg = calllib('mgadgetron', 'mGT_parameter', ...
                self.handle_, 'gadget_chain', id);
            mUtil.checkExecutionStatus(self.name_, hg);
            hv = calllib('mgadgetron', 'mGT_parameter', hg, 'gadget', prop);
            mUtil.checkExecutionStatus(self.name_, hv);
            v = calllib('mutilities', 'mCharDataFromHandle', hv);
            calllib('mutilities', 'mDeleteDataHandle', hg)
            calllib('mutilities', 'mDeleteDataHandle', hv)
        end
    end
end