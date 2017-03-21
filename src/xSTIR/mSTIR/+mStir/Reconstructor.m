classdef Reconstructor < handle
%     Class for generic PET reconstruction objects.
    properties (Constant)
        R = 'Reconstruction';
    end
    properties
        handle
    end
    methods
        function self = Reconstructor()
            self.handle = [];
        end
        function delete(self)
            if ~isempty(self.handle)
                calllib('mstir', 'mDeleteDataHandle', self.handle)
            end
        end
        function set_output_filename_prefix(self, prefix)
            mStir.setParameter(self.handle, self.R, 'output_filename_prefix',...
                prefix, 'c')
        end
    end
end