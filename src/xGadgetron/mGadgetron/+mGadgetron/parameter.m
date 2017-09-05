function value = parameter(handle, set, name, type)
    hv = calllib('mgadgetron', 'mGT_parameter', handle, set, name);
    mUtilities.check_status('parameter', hv)
    if strcmp(type, 'i')
        %value = calllib('mutilities', 'mIntDataFromHandle', hv);
        value = calllib('miutilities', 'mIntDataFromHandle', hv);
    elseif strcmp(type, 'f')
        %value = calllib('mutilities', 'mFloatDataFromHandle', hv);
        value = calllib('miutilities', 'mFloatDataFromHandle', hv);
    end
end