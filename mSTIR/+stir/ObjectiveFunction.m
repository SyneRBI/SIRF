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
            stir.setParameter...
                (self.handle, 'GeneralisedObjectiveFunction', 'prior',...
                prior.handle, 'h')
        end
        function prior = get_prior(self)
            prior = stir.Prior();
            prior.handle = calllib('mstir', 'mSTIR_parameter',...
                self.handle, 'GeneralisedObjectiveFunction', 'prior');
            stir.checkExecutionStatus...
                ('GeneralisedObjectiveFunction:get_prior', prior.handle)
        end
        function set_up(self)
            h = calllib...
                ('mstir', 'mSTIR_setupObject', 'GeneralisedObjectiveFunction',...
                self.handle);
            stir.checkExecutionStatus...
                ('GeneralisedObjectiveFunction:set_up', h)
            calllib('mutilities', 'mDeleteDataHandle', h)
        end
        function v = value(self, image)
            h = calllib('mstir', 'mSTIR_objectiveFunctionValue',...
                self.handle, image.handle);
            stir.checkExecutionStatus...
                ('GeneralisedObjectiveFunction:value', h)
            v = calllib('mutilities', 'mFloatDataFromHandle', h);
            calllib('mutilities', 'mDeleteDataHandle', h)
        end
        function g = gradient(self, image, subset)
            g = stir.Image();
            g.handle = calllib('mstir', 'mSTIR_objectiveFunctionGradient',...
                self.handle, image.handle, subset);
            stir.checkExecutionStatus...
                ('GeneralisedObjectiveFunction:gradient', g.handle)
        end
    end
end