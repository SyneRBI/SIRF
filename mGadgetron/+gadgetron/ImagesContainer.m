classdef ImagesContainer < handle
    properties
        handle_
        name_
    end
    methods
        function self = ImagesContainer(images)
            self.name_ = 'ImagesList';
            self.handle_ = [];
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
            end
        end
        function write(self, file, group)
            handle = calllib('mgadgetron', 'mGT_writeImages', ...
                self.handle_, file, group);
            gadgetron.checkExecutionStatus(self.name_, handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function num = number(self)
            num = calllib('mgadgetron', 'mGT_numImages', self.handle_);
        end
        function r = norm(self)
            handle = calllib('mgadgetron', 'mGT_imagesNorm', ...
                self.handle_);
            gadgetron.checkExecutionStatus('ImagesContainer', handle);
            r = calllib('mgadgetron', 'mDoubleDataFromHandle', handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function z = dot(self, acqs)
            handle = calllib('mgadgetron', 'mGT_imagesDot', ...
                self.handle_, acqs.handle_);
            gadgetron.checkExecutionStatus('ImagesContainer', handle);
            re = calllib('mgadgetron', 'mDoubleReDataFromHandle', handle);
            im = calllib('mgadgetron', 'mDoubleImDataFromHandle', handle);
            z = complex(re, im);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function data = image_as_array(self, im_num)
            ptr_i = libpointer('int32Ptr', zeros(4, 1));
            calllib...
                ('mgadgetron', 'mGT_getImageDimensions', ...
                self.handle_, im_num - 1, ptr_i);
            dim = ptr_i.Value;
            n = dim(1)*dim(2)*dim(3)*dim(4);
            ptr_v = libpointer('doublePtr', zeros(n, 1));
            calllib...
                ('mgadgetron', 'mGT_getImageDataAsDoubleArray', ...
                self.handle_, im_num - 1, ptr_v)
            data = reshape(ptr_v.Value, dim(1), dim(2), dim(3), dim(4));
            %data = reshape(ptr_v.Value, dim(4), dim(3), dim(2), dim(1));
        end
    end
    methods(Static)
        function z = axpby(a, x, b, y)
            z = gadgetron.ImagesContainer();
            z.handle_ = calllib('mgadgetron', 'mGT_imagesZaxpby', ...
                real(a), imag(a), x.handle_, real(b), imag(b), y.handle_);
        end
    end
end