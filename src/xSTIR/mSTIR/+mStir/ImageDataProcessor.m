classdef ImageDataProcessor < handle
    % Class for image filters.
    properties
        name
        handle
    end
    methods
        function self = ImageDataProcessor()
            self.handle = [];
        end
        function delete(self)
            if ~isempty(self.handle)
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
            end
        end
        function apply(self, image)
%         Applies this filter to specified image.
            h = calllib('mstir', 'mSTIR_applyDataProcessor',...
                self.handle, image.handle);
            mUtil.checkExecutionStatus('DataProcessor:apply', h)
            calllib('mutilities', 'mDeleteDataHandle', h)
        end
    end
end