function setParameter(hs, set, par, value, type)
    has_handle = false;
    if ischar(value)
        hv = calllib('mutilities', 'mCharDataHandle', value);
    elseif strcmp(type, 'i')
        hv = calllib('mutilities', 'mIntDataHandle', value);
    elseif strcmp(type, 'f')
        hv = calllib('mutilities', 'mFloatDataHandle', value);
    else
        has_handle = true;
    end
    if has_handle
        handle = calllib...
            ('mstir', 'mSTIR_setParameter', hs, set, par, value.handle);
    else
        handle = calllib('mstir', 'mSTIR_setParameter', hs, set, par, hv);
        calllib('mutilities', 'mDeleteDataHandle', hv)
    end
    mUtil.checkExecutionStatus('setParameter', handle)
    calllib('mutilities', 'mDeleteDataHandle', handle)
end