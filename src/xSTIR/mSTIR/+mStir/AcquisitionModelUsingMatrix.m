classdef AcquisitionModelUsingMatrix < mStir.AcquisitionModel
% ADVANCED USERS ONLY.    
% Class for PET acquisition model with the geometric projection G
% represented by a sparse matrix

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
    end
    methods
        function self = AcquisitionModelUsingMatrix(matrix)
%         Creates an AcquisitionModelUsingMatrix object, optionally setting
%         the ray tracing matrix to be used for projecting;
%         matrix:  a RayTracingMatrix object to represent G in (F).
            self.name = 'AcqModUsingMatrix';
            self.handle = calllib('mstir', 'mSTIR_newObject', self.name);
            mUtil.checkExecutionStatus([self.name ':ctor'], self.handle)
            if nargin < 1
                matrix = mStir.RayTracingMatrix();
            end
            mStir.setParameter...
                (self.handle, self.name, 'matrix', matrix, 'h')
        end
        function delete(self)
            if ~isempty(self.handle)
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
                self.handle = [];
            end
        end
        function set_matrix(self, matrix)
%***SIRF*** set_matrix(matrix) sets the ray tracing matrix to be used for 
%         projecting;
%         matrix:  a RayTracingMatrix object to represent G in (F).
            mStir.setParameter...
                (self.handle, self.name, 'matrix', matrix, 'h')
        end
    end
end