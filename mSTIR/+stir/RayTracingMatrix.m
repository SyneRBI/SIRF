classdef RayTracingMatrix < handle
    properties
        name
        handle
    end
    methods
        function self = RayTracingMatrix()
            self.name = 'RayTracingMatrix';
            self.handle = calllib('mstir', 'mSTIR_newObject', self.name);
        end
        function delete(self)
            calllib('mutilities', 'mDeleteDataHandle', self.handle)
%            calllib('mstir', 'mSTIR_deleteObject', self.handle)
        end
        function set_num_tangential_LORs(self, value)
            stir.setParameter...
                (self.handle, self.name, 'num_tangential_LORs', value, 'i')
        end
        function value = get_num_tangential_LORs(self)
            value = stir.parameter...
                (self.handle, self.name, 'num_tangential_LORs', 'i');
        end
    end
end