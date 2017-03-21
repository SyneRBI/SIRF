classdef AcquisitionInfo < handle
    % Class for acquisition information parameters.
    properties
        flags_
        encode_step_1_
        slice_
        repetition_
    end
    methods
        function self = AcquisitionInfo()
            self.flags_ = 0;
            self.encode_step_1_ = 0;
            self.slice_ = 0;
            self.repetition_ = 0;
        end
    end
end