classdef Gadget < handle
% ADVANCED USERS ONLY. 
% Class for Gadgetron gadgets.

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
    methods (Static)
        function name = class_name()
            name = 'Gadget';
        end
    end
    methods
        function self = Gadget(fullname)
%         Creates a gadget of type and properties specified by the argument,
%         a Matlab string of the form 
%             'gadget_type[(property1=value1[, ...])]'
%         (square brackets embrace optional items, ... stands for etc.).
            [name, prop] = mUtilities.name_and_parameters(fullname);
            self.handle_ = calllib('mgadgetron', 'mGT_newObject', name);
            mUtilities.check_status(name, self.handle_);
            if ~isempty(prop)
                self.set_properties(prop)
            end
            self.name_ = name;
        end
        function delete(self)
            if ~isempty(self.handle_)
                mUtilities.delete(self.handle_)
                %calllib('mutilities', 'mDeleteObject', self.handle_)
            end
        end
        function set_property(self, property, value)
%***SIRF*** set_property(prop, value) assigns value to specified gadget property.
%         prop : property name (Matlab char string)
%         value: property value (Matlab char string)
            handle = calllib('mgadgetron', 'mGT_setGadgetProperty', ...
                self.handle_, property, value);
            mUtilities.check_status(self.name_, handle)
            mUtilities.delete(handle)
            %calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function set_properties(self, properties)
%***SIRF*** Assigns specified values to specified gadget properties.
%         The argument is a Matlab char string with comma-separated list 
%         of property value assignments 'prop_name=prop_value[, ...]'.
            handle = calllib('mgadgetron', 'mGT_setGadgetProperties', ...
                self.handle_, properties);
            mUtilities.check_status(self.name_, handle)
            mUtilities.delete(handle)
            %calllib('mutilities', 'mDeleteDataHandle', handle)
        end
    end
end