classdef RayTracingMatrix < handle
    % Class for objects holding sparse matrix representation of the ray
    % tracing projector G (see AcquisitionModel class).
    properties
        name
        handle
    end
    methods
        function self = RayTracingMatrix()
            self.name = 'RayTracingMatrix';
            self.handle = calllib('mstir', 'mSTIR_newObject', self.name);
            mUtil.checkExecutionStatus(self.name, self.handle)
            mStir.setParameter...
                (self.handle, self.name, 'num_tangential_LORs', 2, 'i')
        end
        function delete(self)
            calllib('mutilities', 'mDeleteDataHandle', self.handle)
        end
        function set_num_tangential_LORs(self, value)
            mStir.setParameter...
                (self.handle, self.name, 'num_tangential_LORs', value, 'i')
        end
        function value = get_num_tangential_LORs(self)
            value = mStir.parameter...
                (self.handle, self.name, 'num_tangential_LORs', 'i');
        end
    end
end