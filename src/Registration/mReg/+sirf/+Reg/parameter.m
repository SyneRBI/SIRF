function value = parameter(handle, set, name, type)
    hv = calllib('mreg', 'mReg_parameter', handle, set, name);
    mUtilities.check_status('parameter', hv)
    if strcmp(type, 'i')
        value = calllib('miutilities', 'mIntDataFromHandle', hv);
    elseif strcmp(type, 'f')
        value = calllib('miutilities', 'mFloatDataFromHandle', hv);
    end
end