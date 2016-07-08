classdef EllipsoidalCylinder < stir.Shape
    properties
        name
    end
    methods
        function self = EllipsoidalCylinder()
            self.name = 'EllipsoidalCylinder';
            self.handle = calllib('mstir', 'mSTIR_newObject', self.name);
        end
        function delete(self)
            if ~isempty(self.handle)
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
%                calllib('mstir', 'mSTIR_deleteObject', self.handle)
                self.handle = [];
            end
        end
        function set_length(self, value)
            stir.setParameter(self.handle, self.name, 'length', value, 'f')
        end
        function value = get_length(self)
            value = stir.parameter(self.handle, self.name, 'length', 'f');
        end
        function set_radii(self, r)
            stir.setParameter(self.handle, self.name, 'radius_x', r(1), 'f')
            stir.setParameter(self.handle, self.name, 'radius_y', r(2), 'f')
        end
    end
end