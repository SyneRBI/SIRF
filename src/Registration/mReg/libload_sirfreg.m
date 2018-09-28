function libload_sirfreg
% load C++-to-C interface library
if ~libisloaded('miutilities')
    fprintf('loading miutilities library...\n')
    try
        [notfound, warnings] = loadlibrary('miutilities');
    catch
        error('miutilities library failed to load\n')
    end
end
% load SIRFReg interface library
if ~libisloaded('msirfreg')
    fprintf('loading msirfreg library...\n')
    try
        [notfound, warnings] = loadlibrary('msirfreg');
    catch
        error('msirfreg library failed to load\n')
    end
end
end