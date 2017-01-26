function value = parameter(handle, set, name, type)
    hv = calllib('mgadgetron', 'mGT_parameter', handle, set, name);
    mGadgetron.checkExecutionStatus('parameter', hv)
    if strcmp(type, 'i')
        value = calllib('mutilities', 'mIntDataFromHandle', hv);
    elseif strcmp(type, 'f')
        value = calllib('mutilities', 'mFloatDataFromHandle', hv);
    end
end