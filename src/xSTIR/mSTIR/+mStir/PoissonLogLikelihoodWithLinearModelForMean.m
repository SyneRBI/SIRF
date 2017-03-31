classdef PoissonLogLikelihoodWithLinearModelForMean < mStir.ObjectiveFunction
%     Class for STIR PoissonLogLikelihoodWithLinearModelForMean object, see
%     http://stir.sourceforge.net/documentation/doxy/html/classstir_1_1PoissonLogLikelihoodWithLinearModelForMean.html
    methods
        function self = PoissonLogLikelihoodWithLinearModelForMean()
            self.handle = [];
        end
        function delete(self)
            if ~isempty(self.handle)
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
                self.handle = [];
            end
        end
        function set_sensitivity_filename(self, name)
            mStir.setParameter...
                (self.handle, 'PoissonLogLikelihoodWithLinearModelForMean',...
                'sensitivity_filename', name, 'c')
        end
        function set_use_subset_sensitivities(self, value)
            if value
                str = 'true';
            else
                str = 'false';
            end
            mStir.setParameter...
                (self.handle, 'PoissonLogLikelihoodWithLinearModelForMean',...
                'use_subset_sensitivities', str, 'c')
        end
        function set_recompute_sensitivity(self, value)
            if value
                str = 'true';
            else
                str = 'false';
            end
            mStir.setParameter...
                (self.handle, 'PoissonLogLikelihoodWithLinearModelForMean',...
                'recompute_sensitivity', str, 'c')
        end
        function sens = get_subset_sensitivity(self, subset)
            sens = mStir.ImageData();
            sens.handle = calllib('mstir', 'mSTIR_subsetSensitivity',...
                self.handle, subset);
            mUtil.checkExecutionStatus...
                ('PoissonLinModMean:get_subset_sensitivity',...
                sens.handle)
        end
        function bar = get_backprojection_of_acquisition_ratio...
                (self, image, subset)
            bar = mStir.ImageData();
            bar.handle = calllib...
                ('mstir', 'mSTIR_objectiveFunctionGradientNotDivided',...
                self.handle, image.handle, subset);
            mUtil.checkExecutionStatus...
                ('PoissonLinModMean:get_backprojection_of_acquisition_ratio',...
                bar.handle)
        end
    end
end
