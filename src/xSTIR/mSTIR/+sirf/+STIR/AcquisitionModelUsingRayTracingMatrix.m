classdef AcquisitionModelUsingRayTracingMatrix < ...
        sirf.STIR.AcquisitionModelUsingMatrix
% Class for PET acquisition model using the ray tracing projection matrix.

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
    end
    methods
        function self = AcquisitionModelUsingRayTracingMatrix(matrix)
%         Creates an AcquisitionModelUsingRayTracingMatrix object.
%         The optional argument sets the ray tracing matrix to be used;
%         matrix:  a RayTracingMatrix object to represent G in (F) -
%                  see AcquisitionModel
            self.name = 'AcqModUsingMatrix';
            self.handle_ = calllib('mstir', 'mSTIR_newObject', self.name);
            sirf.Utilities.check_status([self.name ':ctor'], self.handle_)
            if nargin < 1
                matrix = sirf.STIR.RayTracingMatrix();
            end
            sirf.STIR.setParameter...
                (self.handle_, self.name, 'matrix', matrix, 'h')
        end
        function delete(self)
            if ~isempty(self.handle_)
                %calllib('mutilities', 'mDeleteDataHandle', self.handle_)
                sirf.Utilities.delete(self.handle_)
                self.handle_ = [];
            end
        end
        function matrix = get_matrix(self)
%***SIRF*** Returns the ray tracing matrix used for projecting;
%         matrix:  a RayTracingMatrix object representing G in (F) -
%                  see AcquisitionModel
            matrix = sirf.STIR.RayTracingMatrix();
            sirf.Utilities.delete(matrix.handle_)
            matrix.handle_ = calllib('mstir', 'mSTIR_parameter',...
                self.handle_, self.name, 'matrix');
            sirf.Utilities.check_status...
                ([self.name ':get_matrix'], matrix.handle_)
        end
        function set_num_tangential_LORs(self, value)
%***SIRF*** Set the number of LORs (or rays) for each bin in the sinogram.
%         They are currently (approximately) parallel and spaced in the
%         tangential direction (i.e. orthogonal to the axial direction).
            matrix = self.get_matrix();
            matrix.set_num_tangential_LORs(value)
        end
    end
end