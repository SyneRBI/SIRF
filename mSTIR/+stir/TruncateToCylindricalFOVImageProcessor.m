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
        function set_strictly_less_than_radius(self, flag)
            if flag
                str = 'true';
            else
                str = 'false';
            end
            stir.setParameter(self.handle,...
                'TruncateToCylindricalFOVImageProcessor',...
                'strictly_less_than_radius', str, 'c')
        end
        function flag = get_strictly_less_than_radius(self)
            flag = stir.parameter(self.handle,...
                'TruncateToCylindricalFOVImageProcessor',...
                'strictly_less_than_radius', 'i');
        end
        function apply(self, image)
            h = calllib('mstir', 'mSTIR_applyDataProcessor',...
                self.handle, image.handle);
            stir.checkExecutionStatus...
                ('TruncateToCylindricalFOVImageProcessor:apply', h)
            calllib('mstir', 'mDeleteDataHandle', h)
        end
    end
end