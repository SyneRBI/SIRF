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
                calllib('mstir', 'mDeleteDataHandle', self.handle)
            end
        end
        function apply(self, image)
            h = calllib('mstir', 'mSTIR_applyDataProcessor',...
                self.handle, image.handle);
            stir.checkExecutionStatus('DataProcessor:apply', h)
            calllib('mstir', 'mDeleteDataHandle', h)
        end
    end
end