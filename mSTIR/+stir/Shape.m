classdef Shape < handle
    properties
        handle
    end
    methods
        function self = Shape()
            self.handle = [];
        end
        function delete(self)
            if ~isempty(self.handle)
                calllib('mstir', 'mSTIR_deleteObject', self.handle)
            end
        end
        function set_origin(self, origin)
            stir.setParameter(self.handle, 'Shape', 'x', origin(1), 'f')
            stir.setParameter(self.handle, 'Shape', 'y', origin(2), 'f')
            stir.setParameter(self.handle, 'Shape', 'z', origin(3), 'f')
        end
        function [x, y, z] = get_origin(self)
            x = stir.parameter(self.handle, 'Shape', 'x', 'f');
            y = stir.parameter(self.handle, 'Shape', 'y', 'f');
            z = stir.parameter(self.handle, 'Shape', 'z', 'f');
        end
    end
end