classdef ProjectorsUsingMatrix < handle
    properties
        name
        handle
        matrix
    end
    methods
        function self = ProjectorsUsingMatrix()
            self.name = 'ProjectorsUsingMatrix';
            self.handle = calllib('mstir', 'mSTIR_newObject', self.name);
            self.matrix = [];
        end
        function delete(self)
            calllib('mstir', 'mSTIR_deleteObject', self.handle, 'Projectors')
        end
        function set_matrix(self, matrix)
            stir.setParameter...
                (self.handle, self.name, 'matrix_type', matrix, 'h')
            self.matrix = matrix;
        end
        function matrix = get_matrix(self)
            matrix = self.matrix;
            if isempty(self.matrix)
                error('ProjectorsUsingMatrix:no_matrix_set',...
                    'ProjectorsUsingMatrix: no matrix set')
            end                
        end
    end
end