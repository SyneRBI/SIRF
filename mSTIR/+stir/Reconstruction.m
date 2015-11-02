classdef Reconstruction < handle
    properties
        handle
    end
    methods
        function set_output_filename_prefix(self, prefix)
            stir.setParameter...
                (self.handle, 'Reconstruction', 'output_filename_prefix',...
                prefix, 'c')
        end
    end
end