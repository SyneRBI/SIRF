classdef PoissonLogLh_LinModMean < stir.ObjectiveFunction
    methods
        function self = PoissonLogLh_LinModMean()
            self.handle = [];
        end
        function delete(self)
            if ~isempty(self.handle)
%                calllib('mstir', 'mSTIR_deleteObject', self.handle)
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
                self.handle = [];
            end
        end
        function set_sensitivity_filename(self, name)
            stir.setParameter...
                (self.handle, 'PoissonLogLikelihoodWithLinearModelForMean',...
                'sensitivity_filename', name, 'c')
        end
        function set_use_subset_sensitivities(self, value)
            if value
                str = 'true';
            else
                str = 'false';
            end
            stir.setParameter...
                (self.handle, 'PoissonLogLikelihoodWithLinearModelForMean',...
                'use_subset_sensitivities', str, 'c')
        end
        function set_recompute_sensitivity(self, value)
            if value
                str = 'true';
            else
                str = 'false';
            end
            stir.setParameter...
                (self.handle, 'PoissonLogLikelihoodWithLinearModelForMean',...
                'recompute_sensitivity', str, 'c')
        end
    end
end
