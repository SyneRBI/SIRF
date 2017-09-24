classdef ClientConnector < handle
% ADVANCED USERS ONLY. 
% Class for Gadgetron client connector.

% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC.
% 
% This is software developed for the Collaborative Computational
% Project in Positron Emission Tomography and Magnetic Resonance imaging
% (http://www.ccppetmr.ac.uk/).
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% http://www.apache.org/licenses/LICENSE-2.0
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

    properties
        handle_
        name_
    end
    methods
        function self = ClientConnector()
            self.name_ = 'GTConnector';
            self.handle_ = calllib('mgadgetron', 'mGT_newObject', self.name_);
            mUtilities.check_status(self.name_, self.handle_);
        end
        function delete(self)
            if ~isempty(self.handle_)
                mUtilities.delete(self.handle_)
                %calllib('mutilities', 'mDeleteObject', self.handle_)
            end
        end
        function connect(self, host, port)
            handle = calllib...
                ('mgadgetron', 'mGT_connect', self.handle_, host, port);
            mUtilities.check_status(self.name_, handle);
            mUtilities.delete(handle)
            %calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function disconnect(self)
            handle = calllib...
                ('mgadgetron', 'mGT_disconnect', self.handle_);
            mUtilities.check_status(self.name_, handle);
            mUtilities.delete(handle)
            %calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function register_images_receiver(self, imgs)
            handle = calllib...
                ('mgadgetron', 'mGT_registerImagesReceiver', ...
                self.handle_, imgs.handle_);
            mUtilities.check_status(self.name_, handle);
            mUtilities.delete(handle)
            %calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function config_gadget_chain(self, gc)
            handle = calllib...
                ('mgadgetron', 'mGT_configGadgetChain', ...
                self.handle_, gc.handle_);
            mUtilities.check_status(self.name_, handle);
            mUtilities.delete(handle)
            %calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function send_parameters(self, par)
            handle = calllib...
                ('mgadgetron', 'mGT_sendParameters', self.handle_, par);
            mUtilities.check_status(self.name_, handle);
            mUtilities.delete(handle)
            %calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function send_acquisitions(self, acq)
            handle = calllib...
                ('mgadgetron', 'mGT_sendAcquisitions', ...
                self.handle_, acq.handle_);
            mUtilities.check_status(self.name_, handle);
            mUtilities.delete(handle)
            %calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function send_images(self, img)
            handle = calllib...
                ('mgadgetron', 'mGT_sendImages', self.handle_, img.handle_);
            mUtilities.check_status(self.name_, handle);
            mUtilities.delete(handle)
            %calllib('mutilities', 'mDeleteDataHandle', handle)
        end
    end
end