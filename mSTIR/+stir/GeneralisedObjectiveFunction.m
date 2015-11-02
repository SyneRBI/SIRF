classdef GeneralisedObjectiveFunction < handle
    properties
        name
        handle
        prior
    end
    methods
        function self = GeneralisedObjectiveFunction()
            %fprintf('base class costructor called\n')
            self.handle = calllib('mstir', 'mNewDataHandle');
            self.prior = [];
        end
        function delete(self)
            %fprintf('base class destructor called\n')
            calllib('mstir', 'mDeleteDataHandle', self.handle)
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