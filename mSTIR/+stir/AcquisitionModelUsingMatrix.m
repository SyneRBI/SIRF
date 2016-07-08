classdef AcquisitionModelUsingMatrix < handle
    properties
        name
        handle
        template
        image
    end
    methods
        function self = AcquisitionModelUsingMatrix()
            self.name = 'ProjectorsUsingMatrix';
            self.handle = calllib('mstir', 'mSTIR_newObject', self.name);
            self.template = [];
            self.image = [];
        end
        function delete(self)
            calllib('mutilities', 'mDeleteDataHandle', self.handle)
%            calllib('mstir', 'mSTIR_deleteObject', self.handle)
        end
        function set_matrix(self, matrix)
            stir.setParameter...
                (self.handle, self.name, 'matrix_type', matrix, 'h')
        end
        function matrix = get_matrix(self)
            matrix = stir.RayTracingMatrix();
            matrix.handle = calllib('mstir', 'mSTIR_parameter',...
                self.handle, self.name, 'matrix_type');
            stir.checkExecutionStatus...
                ([self.name ':get_matrix'], matrix.handle)
        end
        function set_up(self, template, image)
            self.template = calllib...
                ('mstir', 'mSTIR_acquisitionModelSetup',...
                self.handle, template, image.handle);
            stir.checkExecutionStatus...
                ([self.name ':set_up'], self.template)
            self.image = calllib('mutilities', 'mCopyOfObjectHandle', image.handle);
%            self.image = calllib('mstir', 'mSTIR_copyOfObject', image.handle);
        end
        function ad = forward(self, image, filename)
            if isempty(self.template)
                error([self.name ':error'],...
                    'forward projection failed: setup not done')
            end
            ad = stir.AcquisitionData();
            ad.handle = calllib('mstir', 'mSTIR_acquisitionModelForward',...
                self.handle, filename, self.template, image.handle);
            stir.checkExecutionStatus...
                ([self.name '.forward'], ad.handle)
            calllib('mutilities', 'mDeleteDataHandle', ad.handle)
            %calllib('mstir', 'mSTIR_deleteObject', ad.handle)
            ad.handle = calllib...
                ('mstir', 'mSTIR_objectFromFile', 'AcquisitionData', filename);
            stir.checkExecutionStatus...
                ([self.name '.forward'], ad.handle)
        end
    end
end