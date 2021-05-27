classdef RayTracingMatrix < handle
% Class for ray tracing matrix objects 
% holding sparse matrix representation of the ray
% tracing projector G (see AcquisitionModel class).

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
        handle_
    end
    methods(Static)
        function name = class_name()
            name = 'RayTracingMatrix';
        end
    end
    methods
        function self = RayTracingMatrix()
            self.name = 'RayTracingMatrix';
            self.handle_ = calllib('mstir', 'mSTIR_newObject', self.name);
            sirf.Utilities.check_status(self.name, self.handle_)
            sirf.STIR.setParameter...
                (self.handle_, self.name, 'num_tangential_LORs', 2, 'i')
        end
        function delete(self)
            sirf.Utilities.delete(self.handle_)
        end
        function set_num_tangential_LORs(self, num)
%***SIRF*** Set the number of LORs (or rays) for each bin in the sinogram.
%         They are currently (approximately) parallel and spaced in the
%         tangential direction (i.e. orthogonal to the axial direction).
            sirf.STIR.setParameter...
                (self.handle_, self.name, 'num_tangential_LORs', num, 'i')
        end
        function value = get_num_tangential_LORs(self)
            value = sirf.STIR.parameter...
                (self.handle_, self.name, 'num_tangential_LORs', 'i');
        end
    end
end