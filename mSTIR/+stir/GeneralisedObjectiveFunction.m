classdef GeneralisedObjectiveFunction < handle
    properties
        name
        handle
        owns_handle
        prior
    end
    methods
        function self = GeneralisedObjectiveFunction()
            self.handle = [];
            self.owns_handle = true;
            self.prior = [];
        end
        function delete(self)
            if self.owns_handle & ~isempty(self.handle)
                calllib('mstir', 'mDeleteDataHandle', self.handle)
            end
        end
        function set_prior(self, prior)
            stir.setParameter...
                (self.handle, 'GeneralisedObjectiveFunction', 'prior',...
                prior.handle, 'h')
            self.prior = prior;
        end
        function prior = get_prior(self)
            prior = self.prior;
            if isempty(self.prior)
                error('GeneralisedObjectiveFunction:no_prior_set',...
                    'GeneralisedObjectiveFunction: no prior set')
            end                
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