classdef GeneralisedPrior < handle
    properties
        name
        handle
    end
    methods
        function self = GeneralisedPrior()
            self.handle = [];
        end
        function delete(self)
            if ~isempty(self.handle)
                calllib('mstir', 'mDeleteDataHandle', self.handle)
            end
        end
        function set_penalisation_factor(self, value)
            stir.setParameter...
                (self.handle, 'GeneralisedPrior', 'penalisation_factor', value, 'f')
        end
        function value = get_penalisation_factor(self)
            value = stir.parameter...
                (self.handle, 'GeneralisedPrior', 'penalisation_factor', 'f');
        end
    end
end