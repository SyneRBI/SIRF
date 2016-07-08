classdef AcquisitionData < handle
    properties
        name
        handle
    end
    methods
        function self = AcquisitionData(filename, template)
            self.handle = [];
            self.name = 'AcquisitionData';
            if nargin < 1
                return
            elseif nargin < 2
                self.handle = calllib...
                    ('mstir', 'mSTIR_objectFromFile',...
                    'AcquisitionData', filename);
            else
                self.handle = calllib...
                    ('mstir', 'mSTIR_acquisitionDataFromTemplate',...
                    filename, template);
            end
            stir.checkExecutionStatus(self.name, self.handle);
        end
        function delete(self)
            if ~isempty(self.handle)
%                calllib('mstir', 'mSTIR_deleteObject', self.handle)
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
            end
        end
        function read_from_file(self, filename)
            if ~isempty(self.handle)
%                calllib('mstir', 'mSTIR_deleteObject', self.handle)
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
            end
            self.handle = calllib...
                ('mstir', 'mSTIR_objectFromFile',...
                'AcquisitionData', filename);
            stir.checkExecutionStatus(self.name, self.handle);
        end
    end
end