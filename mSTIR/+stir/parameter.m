function value = parameter(handle, set, name, type)
    hv = calllib('mstir', 'mSTIR_parameter', handle, set, name);
    stir.checkExecutionStatus('parameter', hv)
    if strcmp(type, 'i')
        value = calllib('mstir', 'mIntDataFromHandle', hv);
    elseif strcmp(type, 'f')
        value = calllib('mstir', 'mFloatDataFromHandle', hv);
    end
end