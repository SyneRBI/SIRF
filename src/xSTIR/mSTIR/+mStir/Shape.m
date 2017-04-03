classdef Shape < handle
% Class for an abstract geometric shape used as a building block for
% creating phantom images.

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

    properties
        handle
    end
    methods
        function self = Shape()
            self.handle = [];
        end
        function delete(self)
            if ~isempty(self.handle)
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
            end
        end
        function set_origin(self, origin)
% ***SIRF*** Sets the (discrete) coordinates of the shape centre on a voxel grid.
            mStir.setParameter(self.handle, 'Shape', 'x', origin(1), 'f')
            mStir.setParameter(self.handle, 'Shape', 'y', origin(2), 'f')
            mStir.setParameter(self.handle, 'Shape', 'z', origin(3), 'f')
        end
        function [x, y, z] = get_origin(self)
% ***SIRF*** Returns the coordinates of the shape centre on a voxel grid.
            x = mStir.parameter(self.handle, 'Shape', 'x', 'f');
            y = mStir.parameter(self.handle, 'Shape', 'y', 'f');
            z = mStir.parameter(self.handle, 'Shape', 'z', 'f');
        end
    end
end