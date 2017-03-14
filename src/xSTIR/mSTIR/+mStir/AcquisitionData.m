classdef AcquisitionData < handle
    % Class for PET acquisition data objects.
    properties
        name
        handle
    end
    methods
        function self = AcquisitionData(arg)
%         Creates new AcquisitionData object from a file or another
%         AcquisitionData object;
%         arg:  file name or AcquisitionData object.
            self.handle = [];
            self.name = 'AcquisitionData';
            if nargin < 1
                return
            elseif ischar(arg)
                self.handle = calllib...
                    ('mstir', 'mSTIR_objectFromFile',...
                    'AcquisitionData', arg);
            elseif isa(arg, 'mStir.AcquisitionData')
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
%         Creates ImageData object containing PET image of dimensions
%         and voxel sizes compatible with the scanner geometry stored
%         in this AcquisitionData object and assigns a given value
%         to all voxels;
%         value: a double.
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
%         Returns 3D array of the acquisition data values.
%         Dimensions are:
%         - number of tangential positions
%         - number of views
%         - number of sinograms
            ptr_i = libpointer('int32Ptr', zeros(3, 1));
            calllib('mstir', 'mSTIR_getAcquisitionsDimensions', ...
                self.handle, ptr_i);
            dim = ptr_i.Value;
            n = dim(1)*dim(2)*dim(3);
            ptr_v = libpointer('doublePtr', zeros(n, 1));
            calllib('mstir', 'mSTIR_getAcquisitionsData', self.handle, ptr_v);
            data = reshape(ptr_v.Value, dim(1), dim(2), dim(3));
        end
        function fill(self, value)
%         Fills the object with values;
%         value: double or array of doubles or an AcquisitionData object
            if isempty(self.handle)
                error([self.name ':fill'], ...
                    'AcquisitionData object not initialized')
            elseif isa(value, 'double')
                if numel(value) > 1
                    ptr_v = libpointer('doublePtr', value);
                    calllib('mstir', 'mSTIR_setAcquisitionsData', ...
                        self.handle, ptr_v);
                else
                    calllib('mstir', 'mSTIR_fillAcquisitionsData', ...
                        self.handle, value);
                end
            elseif isa(value, 'mStir.AcquisitionData')
                calllib('mstir', ...
                    'mSTIR_fillAcquisitionsDataFromAcquisitionsData', ...
                    self.handle, value.handle);
            else
                error([self.name ':fill'], 'wrong fill value')
            end
        end
        function ad = clone(self)
%         Returns a true copy of this object (not Python handle).
            ad = mStir.AcquisitionData(self);
            ad.fill(self)
        end
        function ad = get_empty_copy(self, value)
%         Returns a true copy of this object (not Python handle)
%         filled with a given double value.
            if nargin < 2
                value = 0;
            end
            ad = mStir.AcquisitionData(self);
            ad.fill(value)
        end
    end
end