function value = parameter(handle, set, name, type)
    hv = calllib('mreg', 'mReg_parameter', handle, set, name);
    sirf.Utilities.check_status('parameter', hv)
    if strcmp(type, 'i')
        value = calllib('miutilities', 'mIntDataFromHandle', hv);
    elseif strcmp(type, 'f')
        value = calllib('miutilities', 'mFloatDataFromHandle', hv);
    elseif strcmp(type, 'b')
        value = calllib('miutilities', 'mBoolDataFromHandle', hv);
    end
end