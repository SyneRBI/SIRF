classdef AcquisitionModel < handle
    properties
        handle_
        name_
    end
    methods
        function self = AcquisitionModel(acqs, imgs)
            self.name_ = 'MR_AcquisitionModel';
            self.handle_ = calllib('mgadgetron', 'mGT_AcquisitionModel',...
                acqs.handle_, imgs.handle_);
            mGadgetron.checkExecutionStatus(self.name_, self.handle_);
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
            end
            self.handle_ = [];
        end
        function set_coil_sensitivity_maps(self, csms)
            handle = calllib('mgadgetron', 'mGT_setCSMs', ...
                self.handle_, csms.handle_);
            mGadgetron.checkExecutionStatus(self.name_, handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function acqs = forward(self, imgs)
            acqs = mGadgetron.AcquisitionsContainer();
            acqs.handle_ = calllib...
                ('mgadgetron', 'mGT_AcquisitionModelForward', ...
                self.handle_, imgs.handle_);
            mGadgetron.checkExecutionStatus(self.name_, acqs.handle_);
        end
        function imgs = backward(self, acqs)
            imgs = mGadgetron.ImageData();
            imgs.handle_ = calllib...
                ('mgadgetron', 'mGT_AcquisitionModelBackward', ...
                self.handle_, acqs.handle_);
            mGadgetron.checkExecutionStatus(self.name_, imgs.handle_);
        end
    end
end