classdef Reconstructor < handle
% Class for a generic PET reconstructor object.

% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC.
% Copyright 2015 - 2017 University College London.
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

    properties (Constant)
        R = 'Reconstruction';
    end
    properties
        handle
    end
    methods
        function self = Reconstructor()
            self.handle = [];
        end
        function delete(self)
            if ~isempty(self.handle)
                calllib('mstir', 'mDeleteDataHandle', self.handle)
            end
        end
        function set_output_filename_prefix(self, prefix)
            mStir.setParameter(self.handle, self.R, 'output_filename_prefix',...
                prefix, 'c')
        end
    end
end