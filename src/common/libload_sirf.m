function libload_sirf
% load C++-to-C interface library
if ~libisloaded('miutilities')
    fprintf('loading miutilities library...\n')
    try
        [notfound, warnings] = loadlibrary('miutilities');
    catch
        error('mutilities library failed to load\n')
    end
end
% load STIR interface library
if ~libisloaded('msirf')
    fprintf('loading msirf library...\n')
    try
        [notfound, warnings] = loadlibrary('msirf');
    catch
        error('msirf library failed to load\n')
    end
end
end