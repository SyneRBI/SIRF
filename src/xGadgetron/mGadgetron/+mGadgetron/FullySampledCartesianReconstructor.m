classdef FullySampledCartesianReconstructor < mGadgetron.Reconstructor
% Class to perform reconstruction of fully sampled MR data

    properties
    end
    methods
        function self = FullySampledCartesianReconstructor()
            self.name_ = 'SimpleReconstructionProcessor';
            self.handle_ = calllib('mgadgetron', 'mGT_newObject', self.name_);
            self.input_ = [];
            self.images_ = [];
            mUtil.checkExecutionStatus(self.name_, self.handle_);
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
            end
            self.handle_ = [];
        end
    end
end