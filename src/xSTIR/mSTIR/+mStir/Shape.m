classdef Shape < handle
%     Class for an abstract geometric shape used as a building block for
%     creating phantom images.
    properties
        handle
    end
    methods
        function self = Shape()
            self.handle = [];
        end
        function delete(self)
            if ~isempty(self.handle)
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
            end
        end
        function set_origin(self, origin)
%         Sets the (discrete) coordinates of the shape centre on a voxel grid.
            mStir.setParameter(self.handle, 'Shape', 'x', origin(1), 'f')
            mStir.setParameter(self.handle, 'Shape', 'y', origin(2), 'f')
            mStir.setParameter(self.handle, 'Shape', 'z', origin(3), 'f')
        end
        function [x, y, z] = get_origin(self)
%         Returns the coordinates of the shape centre on a voxel grid.
            x = mStir.parameter(self.handle, 'Shape', 'x', 'f');
            y = mStir.parameter(self.handle, 'Shape', 'y', 'f');
            z = mStir.parameter(self.handle, 'Shape', 'z', 'f');
        end
    end
end