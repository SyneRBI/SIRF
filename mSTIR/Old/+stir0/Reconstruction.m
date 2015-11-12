classdef Reconstruction < handle
    properties (Constant)
        R = 'Reconstruction';
    end
    properties
        handle
    end
    methods
        function self = Reconstruction()
            self.handle = [];
        end
        function delete(self)
            if ~isempty(self.handle)
                calllib('mstir', 'mSTIR_deleteObject', self.handle)
%                calllib('mstir', 'mDeleteDataHandle', self.handle)
            end
        end
        function set_output_filename_prefix(self, prefix)
            stir.setParameter(self.handle, self.R, 'output_filename_prefix',...
                prefix, 'c')
        end
    end
end