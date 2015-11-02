classdef GeneralisedPrior < handle
    properties
        name
        handle
    end
    methods
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