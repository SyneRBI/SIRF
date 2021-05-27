function set_parameter(engine_lib, hs, set, par, value, type)
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
            (engine_lib, 'mSetParameter', hs, set, par, value.handle_);
    else
        handle = calllib(engine_lib, 'mSetParameter', hs, set, par, hv);
        sirf.Utilities.delete(hv)
    end
    sirf.Utilities.check_status('setParameter', handle)
    sirf.Utilities.delete(handle)
end