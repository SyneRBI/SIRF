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
            ('mreg', 'mSetParameter', hs, set, par, value.handle_);
    else
        handle = calllib('mreg', 'mSetParameter', hs, set, par, hv);
        sirf.Utilities.delete(hv)
        %calllib('mutilities', 'mDeleteDataHandle', hv)
    end
    sirf.Utilities.check_status('setParameter', handle)
    sirf.Utilities.delete(handle)
    %calllib('mutilities', 'mDeleteDataHandle', handle)
end