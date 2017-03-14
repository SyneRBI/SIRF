classdef AcquisitionModelUsingMatrix < mStir.AcquisitionModel
    % Class for PET acquisition model with the geometric projection G
    % represented by a sparse matrix
    properties
%         name
%         handle
    end
    methods
        function self = AcquisitionModelUsingMatrix(matrix)
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
            mStir.setParameter...
                (self.handle, self.name, 'matrix', matrix, 'h')
        end
        function matrix = get_matrix(self)
            matrix = mStir.RayTracingMatrix();
            matrix.handle = calllib('mstir', 'mSTIR_parameter',...
                self.handle, self.name, 'matrix');
            mUtil.checkExecutionStatus...
                ([self.name ':get_matrix'], matrix.handle)
        end
    end
end