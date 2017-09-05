classdef GadgetChain < handle
% ADVANCED USERS ONLY. 
% Class for Gadgetron gadget chains.

% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC.
% 
% This is software developed for the Collaborative Computational
% Project in Positron Emission Tomography and Magnetic Resonance imaging
% (http://www.ccppetmr.ac.uk/).
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
            mUtilities.check_status(self.name_, self.handle_);
        end
        function delete(self)
            if ~isempty(self.handle_)
                mUtilities.delete(self.handle_)
                %calllib('mutilities', 'mDeleteObject', self.handle_)
            end
        end
        function add_reader(self, id, reader)
%***SIRF*** add_reader(id, reader) appends the chain with a reader gadget 
%         (a gadget that receives data from the client).
%         id    : gadget id (Matlab string)
%         reader: gadget of reader type (Gadget)
            handle = calllib...
                ('mgadgetron', 'mGT_addReader', self.handle_, id, reader.handle_);
            mUtilities.check_status(self.name_, handle);
            mUtilities.delete(handle)
            %calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function add_writer(self, id, writer)
%***SIRF*** add_writer(id, writer) appends the chain with a writer gadget 
%         (a gadget that sends data to the client).
%         id    : gadget id (Matlab string)
%         writer: gadget of writer type (Gadget)
            handle = calllib...
                ('mgadgetron', 'mGT_addWriter', self.handle_, id, writer.handle_);
            mUtilities.check_status(self.name_, handle);
            mUtilities.delete(handle)
            %calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function add_gadget(self, id, gadget)
%***SIRF*** add_gadget(id, gadget) adds a gadget to the chain.
%         id    : gadget id (Matlab string)
%         gadget: gadget (Gadget)
            handle = calllib...
                ('mgadgetron', 'mGT_addGadget', self.handle_, id, gadget.handle_);
            mUtilities.check_status(self.name_, handle);
            mUtilities.delete(handle)
            %calllib('mutilities', 'mDeleteDataHandle', handle)
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
            mUtilities.check_status(self.name_, hg);
            hv = calllib('mgadgetron', 'mGT_setGadgetProperty', ...
                hg, property, v);
            mUtilities.check_status(self.name_, hv)
            mUtilities.delete(hg)
            mUtilities.delete(hv)
%             calllib('mutilities', 'mDeleteDataHandle', hg)
%             calllib('mutilities', 'mDeleteDataHandle', hv)
        end
        function v = get_gadget_property(self, id, property)
%***SIRF*** get_gadget_property(id, prop) returns the value of the property
%         represented as a Matlab string;
%         id      : gadget id
%         property: property name (Matlab string)
            hg = calllib('mgadgetron', 'mGT_parameter', ...
                self.handle_, 'gadget_chain', id);
            mUtilities.check_status(self.name_, hg);
            hv = calllib('mgadgetron', 'mGT_parameter', hg, 'gadget', property);
            mUtilities.check_status(self.name_, hv);
            %v = calllib('mutilities', 'mCharDataFromHandle', hv);
            v = calllib('miutilities', 'mCharDataFromHandle', hv);
            mUtilities.delete(hg)
            mUtilities.delete(hv)
%             calllib('mutilities', 'mDeleteDataHandle', hg)
%             calllib('mutilities', 'mDeleteDataHandle', hv)
        end
    end
end