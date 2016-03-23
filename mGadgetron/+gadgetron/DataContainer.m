classdef DataContainer < handle
    properties
        handle_
    end
    methods
        function self = DataContainer()
            self.handle_ = [];
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
                self.handle_ = [];
            end
        end
        function num = number(self)
            handle = calllib('mgadgetron', 'mGT_dataItems', self.handle_);
            gadgetron.checkExecutionStatus('AcquisitionsContainer', handle);
            num = calllib('mutilities', 'mIntDataFromHandle', handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function r = norm(self)
            handle = calllib('mgadgetron', 'mGT_norm', self.handle_);
            gadgetron.checkExecutionStatus('AcquisitionsContainer', handle);
            r = calllib('mutilities', 'mDoubleDataFromHandle', handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function z = dot(self, acqs)
            handle = calllib('mgadgetron', 'mGT_dot', self.handle_, ...
                acqs.handle_);
            gadgetron.checkExecutionStatus('AcquisitionsContainer', handle);
            re = calllib('mutilities', 'mDoubleReDataFromHandle', handle);
            im = calllib('mutilities', 'mDoubleImDataFromHandle', handle);
            z = complex(re, im);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
    end
    methods(Static)
        function z = axpby(a, x, b, y)
            z = gadgetron.DataContainer();
            z.handle_ = calllib('mgadgetron', 'mGT_axpby', ...
                real(a), imag(a), x.handle_, real(b), imag(b), y.handle_);
        end
    end
end