classdef DataProcessor < handle
    properties
        name
        handle
    end
    methods
        function self = DataProcessor()
            self.handle = [];
        end
        function delete(self)
            if ~isempty(self.handle)
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
%                calllib('mstir', 'mSTIR_deleteObject', self.handle)
            end
        end
        function apply(self, image)
            h = calllib('mstir', 'mSTIR_applyDataProcessor',...
                self.handle, image.handle);
            stir.checkExecutionStatus('DataProcessor:apply', h)
            calllib('mutilities', 'mDeleteDataHandle', h)
        end
    end
end