classdef ProjectorsUsingMatrix < handle
    properties
        name
        handle
    end
    methods
        function self = ProjectorsUsingMatrix()
            self.name = 'ProjectorsUsingMatrix';
            self.handle = calllib('mstir', 'mSTIR_newObject', self.name);
        end
        function delete(self)
            calllib('mstir', 'mSTIR_deleteObject', self.handle)
        end
        function set_matrix(self, matrix)
            stir.setParameter...
                (self.handle, self.name, 'matrix_type', matrix, 'h')
        end
        function matrix = get_matrix(self)
            matrix = stir.RayTracingMatrix();
            matrix.handle = calllib('mstir', 'mSTIR_parameter',...
                self.handle, self.name, 'matrix_type');
            stir.checkExecutionStatus...
                ([self.name ':get_matrix'], matrix.handle)
        end
    end
end