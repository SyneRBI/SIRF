function libload_reg
% load C++-to-C interface library
if ~libisloaded('miutilities')
    fprintf('loading miutilities library...\n')
    try
        [notfound, warnings] = loadlibrary('miutilities');
    catch
        error('miutilities library failed to load\n')
    end
end
% load Reg interface library
if ~libisloaded('msirfreg')
    fprintf('loading mreg library...\n')
    try
        [notfound, warnings] = loadlibrary('msirfreg');
    catch
        error('mreg library failed to load\n')
    end
end
end