classdef AcquisitionModelUsingMatrix < handle
    properties
        name
        handle
        template
        image
    end
    methods
        function self = AcquisitionModelUsingMatrix()
            self.name = 'AcqModUsingMatrix';
            self.handle = calllib('mstir', 'mSTIR_newObject', self.name);
            stir.checkExecutionStatus([self.name ':ctor'], self.handle)
        end
        function delete(self)
            calllib('mutilities', 'mDeleteDataHandle', self.handle)
        end
        function set_matrix(self, matrix)
            stir.setParameter...
                (self.handle, self.name, 'matrix', matrix, 'h')
        end
        function matrix = get_matrix(self)
            matrix = stir.RayTracingMatrix();
            matrix.handle = calllib('mstir', 'mSTIR_parameter',...
                self.handle, self.name, 'matrix');
            stir.checkExecutionStatus...
                ([self.name ':get_matrix'], matrix.handle)
        end
        function set_up(self, template, image)
            h = calllib...
                ('mstir', 'mSTIR_setupAcquisitionModel',...
                self.handle, template.handle, image.handle);
            stir.checkExecutionStatus([self.name ':set_up'], h)
            %stir.checkExecutionStatus(self.name, h)
            calllib('mutilities', 'mDeleteDataHandle', h)
        end
        function ad = forward(self, image, filename)
            ad = stir.AcquisitionData();
            ad.handle = calllib('mstir', 'mSTIR_acquisitionModelFwd',...
                self.handle, image.handle, filename);
            stir.checkExecutionStatus... %(self.name, ad.handle)
                ([self.name ':forward'], ad.handle)
        end
    end
end