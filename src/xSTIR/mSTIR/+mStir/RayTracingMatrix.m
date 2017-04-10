classdef RayTracingMatrix < handle
% Class for ray tracing matrix objects 
% holding sparse matrix representation of the ray
% tracing projector G (see AcquisitionModel class).

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
        name
        handle
    end
    methods
        function self = RayTracingMatrix()
            self.name = 'RayTracingMatrix';
            self.handle = calllib('mstir', 'mSTIR_newObject', self.name);
            mUtil.checkExecutionStatus(self.name, self.handle)
            mStir.setParameter...
                (self.handle, self.name, 'num_tangential_LORs', 2, 'i')
        end
        function delete(self)
            calllib('mutilities', 'mDeleteDataHandle', self.handle)
        end
        function set_num_tangential_LORs(self, num)
%***SIRF*** Set the number of LORs (or rays) for each bin in the sinogram.
%         They are currently (approximately) parallel and spaced in the
%         tangential direction (i.e. orthogonal to the axial direction).
            mStir.setParameter...
                (self.handle, self.name, 'num_tangential_LORs', num, 'i')
        end
%         function value = get_num_tangential_LORs(self)
%             value = mStir.parameter...
%                 (self.handle, self.name, 'num_tangential_LORs', 'i');
%         end
    end
end