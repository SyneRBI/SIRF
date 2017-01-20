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
%                calllib('mstir', 'mSTIR_deleteObject', self.handle)
                calllib('mutilities', 'mDeleteDataHandle', self.handle)
            end
        end
        function set_prior(self, prior)
            mStir.setParameter...
                (self.handle, 'GeneralisedObjectiveFunction', 'prior',...
                prior.handle, 'h')
        end
        function prior = get_prior(self)
            prior = mStir.Prior();
            prior.handle = calllib('mstir', 'mSTIR_parameter',...
                self.handle, 'GeneralisedObjectiveFunction', 'prior');
            mStir.checkExecutionStatus...
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
            mStir.checkExecutionStatus...
                ('GeneralisedObjectiveFunction:set_up', h)
            calllib('mutilities', 'mDeleteDataHandle', h)
        end
        function v = value(self, image)
            h = calllib('mstir', 'mSTIR_objectiveFunctionValue',...
                self.handle, image.handle);
            mStir.checkExecutionStatus...
                ('GeneralisedObjectiveFunction:value', h)
            v = calllib('mutilities', 'mFloatDataFromHandle', h);
            calllib('mutilities', 'mDeleteDataHandle', h)
        end
        function g = gradient(self, image, subset)
            g = mStir.Image();
            g.handle = calllib('mstir', 'mSTIR_objectiveFunctionGradient',...
                self.handle, image.handle, subset);
            mStir.checkExecutionStatus...
                ('GeneralisedObjectiveFunction:gradient', g.handle)
        end
    end
end