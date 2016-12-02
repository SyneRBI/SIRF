classdef Acquisition < handle
    properties
        handle_
    end
    methods
        function self = Acquisition()
            self.handle_ = [];
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
                self.handle_ = [];
            end
        end
        function ns = number_of_samples(self)
            ns = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'number_of_samples', 'i');
        end
        function ns = active_channels(self)
            ns = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'active_channels', 'i');
        end
    end
end