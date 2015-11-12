classdef TruncateToCylindricalFOVImageProcessor < stir.DataProcessor
    methods
        function self = TruncateToCylindricalFOVImageProcessor(filter)
            self.name = 'TruncateToCylindricalFOVImageProcessor';
            if nargin < 1
                self.handle = calllib...
                    ('mstir', 'mSTIR_newObject', self.name);
            else
                self.handle = calllib...
                    ('mstir', 'mSTIR_copyOfObject', filter.handle);
            end
        end
        function delete(self)
            if ~isempty(self.handle)
                calllib('mstir', 'mSTIR_deleteObject', self.handle)
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