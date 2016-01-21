classdef ClientConnector < handle
    properties
        handle_
        name_
    end
    methods
        function self = ClientConnector()
            self.name_ = 'GTConnector';
            self.handle_ = calllib('mgadgetron', 'mNewObject', self.name_);
            gadgetron.checkExecutionStatus(self.name_, self.handle_);
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mgadgetron', 'mDeleteObject', self.handle_)
            end
        end
        function connect(self, host, port)
            handle = calllib...
                ('mgadgetron', 'mGT_connect', self.handle_, host, port);
            gadgetron.checkExecutionStatus(self.name_, handle);
            calllib('mgadgetron', 'mDeleteDataHandle', handle)
        end
        function disconnect(self)
            handle = calllib...
                ('mgadgetron', 'mGT_disconnect', self.handle_);
            gadgetron.checkExecutionStatus(self.name_, handle);
            calllib('mgadgetron', 'mDeleteDataHandle', handle)
        end
        function register_images_receiver(self, imgs)
            handle = calllib...
                ('mgadgetron', 'mGT_registerImagesReceiver', ...
                self.handle_, imgs.handle_);
            gadgetron.checkExecutionStatus(self.name_, handle);
            calllib('mgadgetron', 'mDeleteDataHandle', handle)
        end
        function config_gadget_chain(self, gc)
            handle = calllib...
                ('mgadgetron', 'mGT_configGadgetChain', ...
                self.handle_, gc.handle_);
            gadgetron.checkExecutionStatus(self.name_, handle);
            calllib('mgadgetron', 'mDeleteDataHandle', handle)
        end
        function send_parameters(self, par)
            handle = calllib...
                ('mgadgetron', 'mGT_sendParameters', self.handle_, par);
            gadgetron.checkExecutionStatus(self.name_, handle);
            calllib('mgadgetron', 'mDeleteDataHandle', handle)
        end
        function send_acquisitions(self, acq)
            handle = calllib...
                ('mgadgetron', 'mGT_sendAcquisitions', ...
                self.handle_, acq.handle_);
            gadgetron.checkExecutionStatus(self.name_, handle);
            calllib('mgadgetron', 'mDeleteDataHandle', handle)
        end
    end
end