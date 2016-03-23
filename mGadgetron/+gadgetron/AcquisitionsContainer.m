classdef AcquisitionsContainer < gadgetron.DataContainer
    properties
    end
    methods
        function self = AcquisitionsContainer()
            self.handle_ = [];
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
                self.handle_ = [];
            end
        end
    end
end