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
                calllib('mstir', 'mSTIR_deleteObject', self.handle)
%                calllib('mstir', 'mDeleteDataHandle', self.handle)
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
            calllib('mstir', 'mDeleteDataHandle', h)
        end
    end
end