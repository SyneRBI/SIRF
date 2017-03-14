classdef AcquisitionModelUsingMatrix < mStir.AcquisitionModel
    % Class for PET acquisition model with the geometric projection G
    % represented by a sparse matrix
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
%         Sets the ray tracing matrix to be used for projecting;
%         matrix:  a RayTracingMatrix object to represent G in (F).
            mStir.setParameter...
                (self.handle, self.name, 'matrix', matrix, 'h')
        end
        function matrix = get_matrix(self)
%         Returns the ray tracing matrix used for projecting;
%         matrix:  a RayTracingMatrix object representing G in (F).
            matrix = mStir.RayTracingMatrix();
            matrix.handle = calllib('mstir', 'mSTIR_parameter',...
                self.handle, self.name, 'matrix');
            mUtil.checkExecutionStatus...
                ([self.name ':get_matrix'], matrix.handle)
        end
    end
end