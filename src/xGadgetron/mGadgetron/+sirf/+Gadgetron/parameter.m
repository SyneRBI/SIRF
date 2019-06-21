function value = parameter(handle, set, name, type, n)
    if nargin < 5
        n = 0;
    end
    value = sirf.Utilities.parameter('mgadgetron', handle, set, name, type, n);
end