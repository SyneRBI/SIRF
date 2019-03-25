classdef AcquisitionModelUsingMatrix < sirf.STIR.AcquisitionModel
% ADVANCED USERS ONLY.    
% Class for PET acquisition model with the geometric projection G
% represented by a sparse matrix - see AcquisitionModel.

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
    end
    methods
        function self = AcquisitionModelUsingMatrix(matrix)
%         Creates an AcquisitionModelUsingMatrix object.
%         The optional argument sets the projection matrix to be used.
%         matrix:  a RayTracingMatrix object to represent G in (F) -
%                  see AcquisitionModel
            self.name = 'AcqModUsingMatrix';
            self.handle_ = calllib('mstir', 'mSTIR_newObject', self.name);
            sirf.Utilities.check_status([self.name ':ctor'], self.handle_)
            if nargin < 1
                matrix = sirf.STIR.RayTracingMatrix();
            end
            sirf.Utilities.assert_validity(matrix, 'RayTracingMatrix')
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
        function set_matrix(self, matrix)
%***SIRF*** set_matrix(matrix) sets the projection matrix to be used.
%         matrix:  a projection matrix object to represent G in (F).
            sirf.Utilities.assert_validity(matrix, 'RayTracingMatrix')
            sirf.STIR.setParameter...
                (self.handle_, self.name, 'matrix', matrix, 'h')
        end
    end
end