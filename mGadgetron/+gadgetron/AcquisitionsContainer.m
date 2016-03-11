classdef AcquisitionsContainer < handle
    properties
        handle_
    end
    methods
        function self = AcquisitionsContainer()
            self.handle_ = [];
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
                self.handle_ = [];
            end
        end
        function r = norm(self)
            handle = calllib('mgadgetron', 'mGT_acquisitionsNorm', ...
                self.handle_);
            gadgetron.checkExecutionStatus('AcquisitionsContainer', handle);
            r = calllib('mgadgetron', 'mDoubleDataFromHandle', handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function z = dot(self, acqs)
            handle = calllib('mgadgetron', 'mGT_acquisitionsDot', ...
                self.handle_, acqs.handle_);
            gadgetron.checkExecutionStatus('AcquisitionsContainer', handle);
            re = calllib('mgadgetron', 'mDoubleReDataFromHandle', handle);
            im = calllib('mgadgetron', 'mDoubleImDataFromHandle', handle);
            z = complex(re, im);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
    end
    methods(Static)
        function z = axpby(a, x, b, y)
            z = gadgetron.AcquisitionsContainer();
            z.handle_ = calllib('mgadgetron', 'mGT_acquisitionsAxpby', ...
                a, x.handle_, b, y.handle_);
        end
    end
end