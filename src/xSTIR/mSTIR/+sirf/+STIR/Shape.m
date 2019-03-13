classdef Shape < handle
% Class for an abstract geometric shape.
% Objects of this class are used as building blocks for creating phantom images.

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
    methods (Static)
        function name = class_name()
            name = 'Shape';
        end
    end
    methods
        function self = Shape()
            self.handle_ = [];
        end
        function delete(self)
            if ~isempty(self.handle_)
                %calllib('mutilities', 'mDeleteDataHandle', self.handle_)
                mUtilities.delete(self.handle_)
            end
        end
        function set_origin(self, origin)
% ***SIRF*** Sets the (discrete) coordinates of the shape centre on a voxel grid.
            sirf.STIR.setParameter(self.handle_, 'Shape', 'x', origin(1), 'f')
            sirf.STIR.setParameter(self.handle_, 'Shape', 'y', origin(2), 'f')
            sirf.STIR.setParameter(self.handle_, 'Shape', 'z', origin(3), 'f')
        end
        function [x, y, z] = get_origin(self)
% ***SIRF*** Returns the coordinates of the shape centre on a voxel grid.
            x = sirf.STIR.parameter(self.handle_, 'Shape', 'x', 'f');
            y = sirf.STIR.parameter(self.handle_, 'Shape', 'y', 'f');
            z = sirf.STIR.parameter(self.handle_, 'Shape', 'z', 'f');
        end
    end
end