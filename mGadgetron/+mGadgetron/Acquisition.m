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
        function f = flags(self)
            f = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'flags', 'i');
        end
        function ns = number_of_samples(self)
            ns = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'number_of_samples', 'i');
        end
        function nc = active_channels(self)
            nc = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'active_channels', 'i');
        end
        function es = idx_kspace_encode_step_1(self)
            es = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'idx_kspace_encode_step_1', 'i');
        end
        function r = idx_repetition(self)
            r = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'idx_repetition', 'i');
        end
        function r = idx_slice(self)
            r = mGadgetron.parameter(self.handle_, ...
                'acquisition', 'idx_slice', 'i');
        end
    end
end