classdef AcquisitionModelData < handle
    properties
        name
        handle
    end
    methods
        function self = AcquisitionModelData(filename)
            self.handle = [];
            self.name = 'AcquisitionModelData';
            if nargin > 0
                self.handle = calllib...
                    ('mstir', 'mSTIR_objectFromFile',...
                    'AcquisitionModelData', filename);
                stir.checkExecutionStatus(self.name, self.handle);
            end
        end
        function delete(self)
            if ~isempty(self.handle)
                calllib('mstir', 'mSTIR_deleteObject', self.handle)
            end
        end
        function read_from_file(self, filename)
            if ~isempty(self.handle)
                calllib('mstir', 'mSTIR_deleteObject', self.handle)
            end
            self.handle = calllib...
                ('mstir', 'mSTIR_objectFromFile',...
                'AcquisitionModelData', filename);
            stir.checkExecutionStatus(self.name, self.handle);
        end
    end
end