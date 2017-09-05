function setParameter(hs, set, par, value, type)
    has_handle = false;
    if ischar(value)
        hv = calllib('miutilities', 'mCharDataHandle', value);
    elseif strcmp(type, 'i')
        hv = calllib('miutilities', 'mIntDataHandle', value);
    elseif strcmp(type, 'f')
        hv = calllib('miutilities', 'mFloatDataHandle', value);
    else
        has_handle = true;
    end
    if has_handle
        handle = calllib...
            ('mstir', 'mSTIR_setParameter', hs, set, par, value.handle);
    else
        handle = calllib('mstir', 'mSTIR_setParameter', hs, set, par, hv);
        mUtilities.delete(hv)
        %calllib('mutilities', 'mDeleteDataHandle', hv)
    end
    mUtilities.check_status('setParameter', handle)
    mUtilities.delete(handle)
    %calllib('mutilities', 'mDeleteDataHandle', handle)
end