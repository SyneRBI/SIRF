classdef Reconstruction < handle
    properties (Constant)
        R = 'Reconstruction';
    end
    properties
        handle
    end
    methods
        function set_output_filename_prefix(self, prefix)
            stir.setParameter(self.handle, self.R, 'output_filename_prefix',...
                prefix, 'c')
        end
    end
end