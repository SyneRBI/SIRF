classdef TruncateToCylindricalFOVImageProcessor < stir.DataProcessor
    properties
        owns_handle
    end
    methods
        function self = TruncateToCylindricalFOVImageProcessor(filter)
            self.name = 'TruncateToCylindricalFOVImageProcessor';
            if nargin < 1
                self.handle = calllib...
                    ('mstir', 'mSTIR_newObject', self.name);
                self.owns_handle = true;
            else
                self.handle = calllib...
                    ('mstir', 'mRefDataHandle', filter.handle);
                self.owns_handle = false;
            end
        end
        function delete(self)
            if self.owns_handle && ~isempty(self.handle)
                calllib('mstir', 'mSTIR_deleteObject', self.handle,...
                    'DataProcessor')
                self.handle = [];
            end
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
    end
end