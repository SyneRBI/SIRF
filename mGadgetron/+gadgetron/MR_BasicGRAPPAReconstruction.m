classdef MR_BasicGRAPPAReconstruction < gadgetron.ImagesReconstructor
    properties
    end
    methods
        function self = MR_BasicGRAPPAReconstruction()
            self.name_ = 'SimpleGRAPPAReconstructionProcessor';
            self.handle_ = calllib('mgadgetron', 'mGT_newObject', self.name_);
            self.input_ = [];
            self.images_ = [];
            gadgetron.checkExecutionStatus(self.name_, self.handle_);
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
        function [images, gfactors] = get_output(self)
            gf = self.value_of_gadget_property('gadget4', 'send_out_gfactor');
            output = get_output@gadgetron.ImagesReconstructor(self);
            if strcmp(gf, 'true')
                images = output.select(2);
                gfactors = output.select(2, 1);
            else
                images = output;
                gfactors = [];
            end
        end
    end
end