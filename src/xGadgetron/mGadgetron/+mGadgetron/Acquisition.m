classdef Acquisition < handle
% INTERNAL USE ONLY.
% Class for the ISMRMRD acquisition object.

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
    end
    methods
        function self = Acquisition()
            self.handle_ = [];
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
                self.handle_ = [];
            end
        end
        function f = flags(self)
            f = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'flags', 'i');
        end
        function ns = number_of_samples(self)
            ns = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'number_of_samples', 'i');
        end
        function nc = active_channels(self)
            nc = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'active_channels', 'i');
        end
        function es = idx_kspace_encode_step_1(self)
            es = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'idx_kspace_encode_step_1', 'i');
        end
        function r = idx_repetition(self)
            r = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'idx_repetition', 'i');
        end
        function r = idx_slice(self)
            r = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'idx_slice', 'i');
        end
    end
end