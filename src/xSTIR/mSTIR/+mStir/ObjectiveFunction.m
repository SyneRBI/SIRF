classdef ObjectiveFunction < handle
    properties
        name
        handle
    end
    methods
        function self = ObjectiveFunction()
            self.handle = [];
        end
        function delete(self)
            if ~isempty(self.handle)
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
                self.handle = [];
            end
        end
        function set_prior(self, prior)
            mStir.setParameter...
                (self.handle, 'GeneralisedObjectiveFunction', 'prior',...
                prior, 'h')
        end
        function prior = get_prior(self)
            prior = mStir.Prior();
            prior.handle = calllib('mstir', 'mSTIR_parameter',...
                self.handle, 'GeneralisedObjectiveFunction', 'prior');
            mUtil.checkExecutionStatus...
                ('GeneralisedObjectiveFunction:get_prior', prior.handle)
        end
        function set_num_subsets(self, num)
            mStir.setParameter...
                (self.handle, 'GeneralisedObjectiveFunction', ...
                'num_subsets', num, 'i')
        end
        function set_up(self, image)
            h = calllib('mstir', 'mSTIR_setupObjectiveFunction', ...
                self.handle, image.handle);
%                 ('mstir', 'mSTIR_setupObject', 'GeneralisedObjectiveFunction',...
%                 self.handle);
            mUtil.checkExecutionStatus...
                ('GeneralisedObjectiveFunction:set_up', h)
            calllib('mutilities', 'mDeleteDataHandle', h)
        end
        function v = value(self, image)
            h = calllib('mstir', 'mSTIR_objectiveFunctionValue',...
                self.handle, image.handle);
            mUtil.checkExecutionStatus...
                ('GeneralisedObjectiveFunction:value', h)
            v = calllib('mutilities', 'mFloatDataFromHandle', h);
            calllib('mutilities', 'mDeleteDataHandle', h)
        end
        function g = gradient(self, image, subset)
            g = mStir.ImageData();
            g.handle = calllib('mstir', 'mSTIR_objectiveFunctionGradient',...
                self.handle, image.handle, subset);
            mUtil.checkExecutionStatus...
                ('GeneralisedObjectiveFunction:gradient', g.handle)
        end
    end
end