classdef OSMAPOSLReconstruction < stir.IterativeReconstruction
    properties
        name
    end
    methods
        function self = OSMAPOSLReconstruction(filename)
            self.name = 'OSMAPOSL';
            if nargin > 0
                self.handle = calllib('mstir', 'mSTIR_newReconstruction',...
                    'OSMAPOSL', filename);
            else
                self.handle = calllib('mstir', 'mSTIR_newReconstruction',...
                    'OSMAPOSL', '');
            end
            stir.checkExecutionStatus(self.name, self.handle);
        end
        function delete(self)
            calllib('mstir', 'mSTIR_deleteReconstruction', self.handle)
        end
        function set_MAP_model(self, model)
            stir.setParameter(self.handle, self.name, 'MAP_model', model, 'c')
        end
    end
end
