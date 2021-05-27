function value = parameter(engine_lib, handle, set, name, type, n)
    hv = calllib(engine_lib, 'mParameter', handle, set, name);
    sirf.Utilities.check_status('parameter', hv)
    if strcmp(type, 'i')
        value = calllib('miutilities', 'mIntDataFromHandle', hv);
    elseif strcmp(type, 'f')
        value = calllib('miutilities', 'mFloatDataFromHandle', hv);
    elseif strcmp(type, 'b')
        value = calllib('miutilities', 'mBoolDataFromHandle', hv);
    elseif strcmp(type, 'u16s')
        value = zeros(n, 1);
        for i = 1 : n
            value(i) = calllib('miutilities', 'mUint16DataItemFromHandle', ...
                hv, i - 1);
        end
    elseif strcmp(type, 'u32s')
        value = zeros(n, 1);
        for i = 1 : n
            value(i) = calllib('miutilities', 'mUint32DataItemFromHandle', ...
                hv, i - 1);
        end
    elseif strcmp(type, 'u64s')
        value = zeros(n, 1);
        for i = 1 : n
            value(i) = calllib('miutilities', 'mUint64DataItemFromHandle', ...
                hv, i - 1);
        end
    elseif strcmp(type, 'fs')
        value = zeros(n, 1);
        for i = 1 : n
            value(i) = calllib('miutilities', 'mFloatDataItemFromHandle', ...
                hv, i - 1);
        end
    end
end