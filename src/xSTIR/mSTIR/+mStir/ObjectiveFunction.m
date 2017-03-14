classdef ObjectiveFunction < handle
    % Class for the objective function maximized by the iterative reconstruction
    % algorithms.
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
%         Sets the prior (penalty term to be added to the objective function).
            mStir.setParameter...
                (self.handle, 'GeneralisedObjectiveFunction', 'prior',...
                prior, 'h')
        end
        function prior = get_prior(self)
%         Returns the prior currently used by this objective function.
            prior = mStir.Prior();
            prior.handle = calllib('mstir', 'mSTIR_parameter',...
                self.handle, 'GeneralisedObjectiveFunction', 'prior');
            mUtil.checkExecutionStatus...
                ('GeneralisedObjectiveFunction:get_prior', prior.handle)
        end
        function set_num_subsets(self, num)
%         Sets the number of subsets of ray projections to be used for computing
%         additive components of the gradient used by Ordered Subset algorithms
%         for maximizing this objective function.
%         Assuming for simplicity of illustration that the ray tracing projector G
%         is a matrix, the subsets in question are subsets of its rows.
            mStir.setParameter...
                (self.handle, 'GeneralisedObjectiveFunction', ...
                'num_subsets', num, 'i')
        end
        function set_up(self, image)
%         Prepares this object for use.
            h = calllib('mstir', 'mSTIR_setupObjectiveFunction', ...
                self.handle, image.handle);
            mUtil.checkExecutionStatus('GeneralisedObjectiveFunction:set_up', h)
            calllib('mutilities', 'mDeleteDataHandle', h)
        end
        function v = value(self, image)
%         Returns the value of this objective function on <image>.
            h = calllib('mstir', 'mSTIR_objectiveFunctionValue',...
                self.handle, image.handle);
            mUtil.checkExecutionStatus...
                ('GeneralisedObjectiveFunction:value', h)
            v = calllib('mutilities', 'mFloatDataFromHandle', h);
            calllib('mutilities', 'mDeleteDataHandle', h)
        end
        function g = gradient(self, image, subset)
%         Returns the value of the additive component of the gradient of this 
%         objective function on <image> corresponding to a given <subset>
%         (see set_num_subsets() method).
            g = mStir.ImageData();
            g.handle = calllib('mstir', 'mSTIR_objectiveFunctionGradient',...
                self.handle, image.handle, subset);
            mUtil.checkExecutionStatus...
                ('GeneralisedObjectiveFunction:gradient', g.handle)
        end
    end
end