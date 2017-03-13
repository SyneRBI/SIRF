classdef AcquisitionData < handle
    % Class for MR acquisition data objects.
    properties
        name
        handle
    end
    methods
        function self = AcquisitionData(arg)
            self.handle = [];
            self.name = 'AcquisitionData';
            if nargin < 1
                return
            elseif ischar(arg)
                self.handle = calllib...
                    ('mstir', 'mSTIR_objectFromFile',...
                    'AcquisitionData', arg);
            elseif isa(arg, mStir.AcquisitionData)
                self.handle = calllib...
                    ('mstir', 'mSTIR_acquisitionsDataFromTemplate',...
                    arg.handle);
            else
                error('AcquisitionData:wrong_ctor_source', ...
                'wrong source in AcquisitionData constructor')
            end
            mUtil.checkExecutionStatus(self.name, self.handle);
        end
        function delete(self)
            if ~isempty(self.handle)
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
            end
        end
        function read_from_file(self, filename)
            if ~isempty(self.handle)
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
            end
            self.handle = calllib('mstir', 'mSTIR_objectFromFile', ...
                'AcquisitionData', filename);
            mUtil.checkExecutionStatus(self.name, self.handle);
        end
        function image = create_empty_image(self, value)
            image = mStir.ImageData();
            image.handle = calllib...
                ('mstir', 'mSTIR_imageFromAcquisitionData', self.handle);
            mUtil.checkExecutionStatus...
                ([self.name ':create_empty_image'], image.handle);
            if nargin > 1
                image.fill(value)
            end
        end
        function data = as_array(self)
%         Returns 3D array of this acquisition data values.

            ptr_i = libpointer('int32Ptr', zeros(3, 1));
            calllib('mstir', 'mSTIR_getAcquisitionsDimensions', ...
                self.handle, ptr_i);
            dim = ptr_i.Value;
            n = dim(1)*dim(2)*dim(3);
            ptr_v = libpointer('doublePtr', zeros(n, 1));
            calllib('mstir', 'mSTIR_getAcquisitionsData', self.handle, ptr_v);
            data = reshape(ptr_v.Value, dim(1), dim(2), dim(3));
        end
    end
end