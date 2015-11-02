classdef TruncateToCylindricalFOVImageProcessor < handle
    properties
        name
        handle
    end
    methods
        function self = TruncateToCylindricalFOVImageProcessor()
            self.name = 'TruncateToCylindricalFOVImageProcessor';
            self.handle = calllib('mstir', 'mSTIR_newObject', self.name);
        end
        function delete(self)
            calllib('mstir', 'mSTIR_deleteObject', self.handle, 'DataProcessor')
        end
    end
end