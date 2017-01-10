classdef ImagesContainer < mGadgetron.DataContainer
    properties
        name_
    end
    methods
        function self = ImagesContainer()
            self.name_ = 'ImagesContainer';
            self.handle_ = [];
        end
        function delete(self)
            if ~isempty(self.handle_)
                calllib('mutilities', 'mDeleteObject', self.handle_)
                self.handle_ = [];
            end
        end
        function write(self, file, group)
            handle = calllib('mgadgetron', 'mGT_writeImages', ...
                self.handle_, file, group);
            mGadgetron.checkExecutionStatus(self.name_, handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
        end
        function images = select(self, attr, value)
            images = mGadgetron.ImagesContainer();
            images.handle_ = calllib('mgadgetron', 'mGT_selectImages', ...
                self.handle_, attr, value);
            mGadgetron.checkExecutionStatus(self.name_, images.handle_);
        end
        function images = process(self, list)
            ip = mGadgetron.ImagesProcessor(list);
            images = ip.process(self);
        end
        function ft = is_real(self)
            handle = calllib('mgadgetron', 'mGT_imageDataType', ...
                self.handle_, 0);
            mGadgetron.checkExecutionStatus(self.name_, handle);
            v = calllib('mutilities', 'mIntDataFromHandle', handle);
            calllib('mutilities', 'mDeleteDataHandle', handle)
            ft = (v ~= 7 && v ~= 8);
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
        end
        function data = as_array(self)
            ptr_i = libpointer('int32Ptr', zeros(4, 1));
            calllib('mgadgetron', 'mGT_getImageDimensions', ...
                self.handle_, 0, ptr_i);
            dim = ptr_i.Value;
            nz = dim(3)*dim(4)*self.number();
            n = dim(1)*dim(2)*nz;
            if self.is_real()
                ptr_v = libpointer('doublePtr', zeros(n, 1));
                calllib('mgadgetron', 'mGT_getImagesDataAsDoubleArray', ...
                    self.handle_, ptr_v)
                data = reshape(ptr_v.Value, dim(1), dim(2), nz);
            else
                ptr_re = libpointer('doublePtr', zeros(n, 1));
                ptr_im = libpointer('doublePtr', zeros(n, 1));
                calllib...
                    ('mgadgetron', 'mGT_getImagesDataAsComplexArray', ...
                    self.handle_, ptr_re, ptr_im)
                re = reshape(ptr_re.Value, dim(1), dim(2), nz);
                im = reshape(ptr_im.Value, dim(1), dim(2), nz);
                data = re + 1j*im;
            end
        end
        function show(self)
            ni = self.number();
            if ni < 1
                return
            end
            data = self.as_array();
            if ~self.is_real()
                data = abs(data);
            end
            data = data/max(max(max(data)));
            fprintf('Please enter the number of the image to view\n')
            fprintf('(a value outside the range [1 : %d] will stop this loop)\n', ni)
            while true
                i = input('image: ');
                if i < 1 || i > ni
                    break
                end
                figure(i)
                imshow(data(:, :, i))
                title(['image ' num2str(i)])
            end
        end
    end
end