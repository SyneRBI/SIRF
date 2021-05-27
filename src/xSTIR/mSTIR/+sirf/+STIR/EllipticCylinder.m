classdef EllipticCylinder < sirf.STIR.Shape
% Class for elliptic cylinder shape.

% SyneRBI Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC.
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
        name
    end
    methods
        function self = EllipticCylinder()
%         Creates an EllipticCylinder object.
            self.name = 'EllipsoidalCylinder';
            self.handle_ = calllib('mstir', 'mSTIR_newObject', self.name);
        end
        function delete(self)
            if ~isempty(self.handle_)
                %calllib('mutilities', 'mDeleteDataHandle', self.handle_)
                sirf.Utilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function set_length(self, value)
%***SIRF*** Sets the length (height) of the cylinder.
            sirf.STIR.setParameter(self.handle_, self.name, 'length', value, 'f')
        end
        function value = get_length(self)
%***SIRF*** Returns the length (height) of the cylinder.
            value = sirf.STIR.parameter(self.handle_, self.name, 'length', 'f');
        end
        function set_radii(self, r)
%***SIRF*** Sets the radii of the cylinder.
            sirf.STIR.setParameter(self.handle_, self.name, 'radius_x', r(1), 'f')
            sirf.STIR.setParameter(self.handle_, self.name, 'radius_y', r(2), 'f')
        end
    end
end