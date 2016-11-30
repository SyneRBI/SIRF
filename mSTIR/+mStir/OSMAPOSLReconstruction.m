classdef OSMAPOSLReconstruction < mStir.IterativeReconstruction
    properties
        name
    end
    methods
        function self = OSMAPOSLReconstruction(filename)
            self.name = 'OSMAPOSL';
            if nargin < 1
                filename = '';
            end
            self.handle = calllib...
                ('mstir', 'mSTIR_objectFromFile',...
                'OSMAPOSLReconstruction', filename);
            mStir.checkExecutionStatus(self.name, self.handle);
        end
        function delete(self)
            calllib('mutilities', 'mDeleteDataHandle', self.handle)
%            calllib('mstir', 'mSTIR_deleteObject', self.handle)
            self.handle = [];
        end
        function set_MAP_model(self, model)
            mStir.setParameter(self.handle, self.name, 'MAP_model', model, 'c')
        end
        function obj_fun = get_objective_function(self)
            obj_fun = mStir.PoissonLogLh_LinModMean();
            obj_fun.handle = calllib('mstir', 'mSTIR_parameter',...
                self.handle, self.name, 'objective_function');
            mStir.checkExecutionStatus...
                ([self.name ':get_objective_function'], obj_fun.handle)
        end
    end
end
