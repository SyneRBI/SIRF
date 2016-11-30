classdef ClientConnector < handle
    properties
        handle_
        name_
    end
    methods
        function self = ClientConnector()
            self.name_ = 'GTConnector';
            self.handle_ = calllib('mgadgetron', 'mGT_newObject', self.name_);
            mGadgetron.checkExecutionStatus(self.name_, self.handle_);
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
            end
        end
        function connect(self, host, port)
            handle = calllib...
                ('mgadgetron', 'mGT_connect', self.handle_, host, port);
            mGadgetron.checkExecutionStatus(self.name_, handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function disconnect(self)
            handle = calllib...
                ('mgadgetron', 'mGT_disconnect', self.handle_);
            mGadgetron.checkExecutionStatus(self.name_, handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function register_images_receiver(self, imgs)
            handle = calllib...
                ('mgadgetron', 'mGT_registerImagesReceiver', ...
                self.handle_, imgs.handle_);
            mGadgetron.checkExecutionStatus(self.name_, handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function config_gadget_chain(self, gc)
            handle = calllib...
                ('mgadgetron', 'mGT_configGadgetChain', ...
                self.handle_, gc.handle_);
            mGadgetron.checkExecutionStatus(self.name_, handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function send_parameters(self, par)
            handle = calllib...
                ('mgadgetron', 'mGT_sendParameters', self.handle_, par);
            mGadgetron.checkExecutionStatus(self.name_, handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function send_acquisitions(self, acq)
            handle = calllib...
                ('mgadgetron', 'mGT_sendAcquisitions', ...
                self.handle_, acq.handle_);
            mGadgetron.checkExecutionStatus(self.name_, handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function send_images(self, img)
            handle = calllib...
                ('mgadgetron', 'mGT_sendImages', self.handle_, img.handle_);
            mGadgetron.checkExecutionStatus(self.name_, handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
    end
end