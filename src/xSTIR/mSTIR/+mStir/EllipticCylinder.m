classdef EllipticCylinder < mStir.Shape
% Class for ellipsoidal cylinder shape.
    properties
        name
    end
    methods
        function self = EllipticCylinder()
            self.name = 'EllipsoidalCylinder';
            self.handle = calllib('mstir', 'mSTIR_newObject', self.name);
        end
        function delete(self)
            if ~isempty(self.handle)
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
                self.handle = [];
            end
        end
        function set_length(self, value)
            mStir.setParameter(self.handle, self.name, 'length', value, 'f')
        end
        function value = get_length(self)
            value = mStir.parameter(self.handle, self.name, 'length', 'f');
        end
        function set_radii(self, r)
            mStir.setParameter(self.handle, self.name, 'radius_x', r(1), 'f')
            mStir.setParameter(self.handle, self.name, 'radius_y', r(2), 'f')
        end
    end
end