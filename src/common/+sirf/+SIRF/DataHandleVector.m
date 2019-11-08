classdef DataHandleVector < handle
% INTERNAL USE ONLY.
% Class for an abstract data container.

% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2019 University College London
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
        name
        handle_
    end
    methods(Static)
        function name = class_name()
            name = 'DataHandleVector';
        end
    end
    methods
        function self = DataHandleVector()
            self.name = 'DataHandleVector';
            self.handle_ = calllib('msirf', 'mSIRF_newObject', self.name);
            sirf.Utilities.check_status(self.name, self.handle_);
        end
        function delete(self)
            if ~isempty(self.handle_)
                sirf.Utilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function push_back(self,handle_)
        	%Push back new data handle.
        	calllib('msirf', 'mSIRF_DataHandleVector_push_back', self.handle_, handle_);
        end
    end
end