classdef ImageFilter < handle
    % Abstract base class for image filters.
    properties
        name
        handle
    end
    methods
        function self = ImageFilter()
            self.handle = [];
        end
        function delete(self)
            if ~isempty(self.handle)
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
            end
        end
        function apply(self, image)
%         Applies this filter to <image>.
            h = calllib('mstir', 'mSTIR_applyDataProcessor',...
                self.handle, image.handle);
            mUtil.checkExecutionStatus('DataProcessor:apply', h)
            calllib('mutilities', 'mDeleteDataHandle', h)
        end
    end
end