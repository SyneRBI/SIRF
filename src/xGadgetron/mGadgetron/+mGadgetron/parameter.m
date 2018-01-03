function value = parameter(handle, set, name, type, n)
    hv = calllib('mgadgetron', 'mGT_parameter', handle, set, name);
    mUtilities.check_status('parameter', hv)
    if strcmp(type, 'i')
        value = calllib('miutilities', 'mIntDataFromHandle', hv);
    elseif strcmp(type, 'f')
        value = calllib('miutilities', 'mFloatDataFromHandle', hv);
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