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
        function obj = same_object(self)
            if isa(self, class(mGadgetron.ImageData))
                obj = mGadgetron.ImageData();
            elseif isa(self, class(mGadgetron.AcquisitionsContainer))
                obj = mGadgetron.AcquisitionsContainer();
            else
                obj = mGadgetron.DataContainer();
            end
        end
        function num = number(self)
            handle = calllib('mgadgetron', 'mGT_dataItems', self.handle_);
            mGadgetron.checkExecutionStatus('DataContainer', handle);
            num = calllib('mutilities', 'mIntDataFromHandle', handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function r = norm(self)
            handle = calllib('mgadgetron', 'mGT_norm', self.handle_);
            mGadgetron.checkExecutionStatus('DataContainer', handle);
            r = calllib('mutilities', 'mDoubleDataFromHandle', handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function z = dot(self, other)
            handle = calllib('mgadgetron', 'mGT_dot', self.handle_, ...
                other.handle_);
            mGadgetron.checkExecutionStatus('DataContainer', handle);
            re = calllib('mutilities', 'mDoubleReDataFromHandle', handle);
            im = calllib('mutilities', 'mDoubleImDataFromHandle', handle);
            z = complex(re, im);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function z = minus(self, other)
            z = self.same_object();
            z.handle_ = calllib('mgadgetron', 'mGT_axpby', ...
                1.0, 0.0, self.handle_, -1.0, 0.0, other.handle_);
        end
        function z = mtimes(self, other)
            if isobject(other)
                z = self.dot(other);
            elseif isreal(other)
                z = self.same_object();
                z.handle_ = calllib('mgadgetron', 'mGT_axpby', ...
                    other, 0.0, self.handle_, 0.0, 0.0, self.handle_);
            else
                z = self.same_object();
                z.handle_ = calllib('mgadgetron', 'mGT_axpby', ...
                    real(other), imag(other), self.handle_, ...
                    0.0, 0.0, self.handle_);
            end
        end
    end
    methods(Static)
        function z = axpby(a, x, b, y)
            z = self.same_object();
            z.handle_ = calllib('mgadgetron', 'mGT_axpby', ...
                real(a), imag(a), x.handle_, real(b), imag(b), y.handle_);
        end
    end
end