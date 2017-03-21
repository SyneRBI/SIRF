classdef CartesianGRAPPAReconstructor < mGadgetron.Reconstructor
    % Class for a reconstructor from undersampled Cartesian raw data.
    properties
    end
    methods
        function self = CartesianGRAPPAReconstructor()
            self.name_ = 'SimpleGRAPPAReconstructionProcessor';
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
        function compute_gfactors(self, flag)
            self.set_gadget_property('gadget4', 'send_out_gfactor', flag)
        end
    end
end