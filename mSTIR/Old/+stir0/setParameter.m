function setParameter(hs, set, par, value, type)
    has_handle = false;
    if ischar(value)
        hv = calllib('mstir', 'mCharDataHandle', value);
    elseif strcmp(type, 'i')
        hv = calllib('mstir', 'mIntDataHandle', value);
    elseif strcmp(type, 'f')
        hv = calllib('mstir', 'mFloatDataHandle', value);
    else
        has_handle = true;
    end
    if has_handle
        handle = calllib...
            ('mstir', 'mSTIR_setParameter', hs, set, par, value.handle);
    else
        handle = calllib('mstir', 'mSTIR_setParameter', hs, set, par, hv);
        calllib('mstir', 'mDeleteDataHandle', hv)
    end
    stir.checkExecutionStatus('setParameter', handle)
    calllib('mstir', 'mDeleteDataHandle', handle)
end