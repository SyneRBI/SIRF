classdef GadgetChain < handle
% ADVANCED USERS ONLY. 
% Class for Gadgetron gadget chains.

% SyneRBI Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2020 Rutherford Appleton Laboratory STFC.
% 
% This is software developed for the Collaborative Computational
% Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
% (http://www.ccpsynerbi.ac.uk/).
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% http://www.apache.org/licenses/LICENSE-2.0
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

    properties
        handle_
        name_
    end
    methods
        function self = GadgetChain()
%         Creates an empty Gadgetron chain.
            self.name_ = 'GadgetChain';
            self.handle_ = calllib('mgadgetron', 'mGT_newObject', self.name_);
            sirf.Utilities.check_status(self.name_, self.handle_);
        end
        function delete(self)
            if ~isempty(self.handle_)
                sirf.Utilities.delete(self.handle_)
            end
        end
%         function add_reader(self, id, reader)
% %***SIRF*** add_reader(id, reader) appends the chain with a reader gadget 
% %         (a gadget that receives data from the client).
% %         id    : gadget id (Matlab string)
% %         reader: gadget of reader type (Gadget)
%             handle = calllib...
%                 ('mgadgetron', 'mGT_addReader', self.handle_, id, reader.handle_);
%             sirf.Utilities.check_status(self.name_, handle);
%             sirf.Utilities.delete(handle)
%         end
%         function add_writer(self, id, writer)
% %***SIRF*** add_writer(id, writer) appends the chain with a writer gadget 
% %         (a gadget that sends data to the client).
% %         id    : gadget id (Matlab string)
% %         writer: gadget of writer type (Gadget)
%             handle = calllib...
%                 ('mgadgetron', 'mGT_addWriter', self.handle_, id, writer.handle_);
%             sirf.Utilities.check_status(self.name_, handle);
%             sirf.Utilities.delete(handle)
%         end
        function set_host(self, host)
            handle = calllib('mgadgetron', 'mGT_setHost', self.handle_, host);
            sirf.Utilities.check_status(self.name_, handle);
            sirf.Utilities.delete(handle)
        end
        function set_port(self, port)
            handle = calllib('mgadgetron', 'mGT_setPort', self.handle_, port);
            sirf.Utilities.check_status(self.name_, handle);
            sirf.Utilities.delete(handle)
        end
        function add_gadget(self, id, gadget)
%***SIRF*** add_gadget(id, gadget) adds a gadget to the chain.
%         id    : gadget id (Matlab string)
%         gadget: gadget (Gadget)
            sirf.Utilities.assert_validity(gadget, 'Gadget')
            handle = calllib...
                ('mgadgetron', 'mGT_addGadget', self.handle_, id, gadget.handle_);
            sirf.Utilities.check_status(self.name_, handle);
            sirf.Utilities.delete(handle)
        end
        function set_gadget_property(self, id, property, value)
%***SIRF*** set_gadget_property(id, prop, val) assigns value to gadget property.
%         id   : gadget id
%         prop : property name (Matlab string)
%         value: property value
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
            sirf.Utilities.check_status(self.name_, hg);
            hv = calllib('mgadgetron', 'mGT_setGadgetProperty', ...
                hg, property, v);
            sirf.Utilities.check_status(self.name_, hv)
            sirf.Utilities.delete(hg)
            sirf.Utilities.delete(hv)
        end
        function v = get_gadget_property(self, id, property)
%***SIRF*** get_gadget_property(id, prop) returns the value of the property
%         represented as a Matlab string;
%         id      : gadget id
%         property: property name (Matlab string)
            hg = calllib('mgadgetron', 'mGT_parameter', ...
                self.handle_, 'gadget_chain', id);
            sirf.Utilities.check_status(self.name_, hg);
            hv = calllib('mgadgetron', 'mGT_parameter', hg, 'gadget', property);
            sirf.Utilities.check_status(self.name_, hv);
            %v = calllib('mutilities', 'mCharDataFromHandle', hv);
            v = calllib('miutilities', 'mCharDataFromHandle', hv);
            sirf.Utilities.delete(hg)
            sirf.Utilities.delete(hv)
        end
    end
end
