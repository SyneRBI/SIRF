function libload_stir
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
if ~libisloaded('mstir')
    fprintf('loading mstir library...\n')
    try
        [notfound, warnings] = loadlibrary('mstir');
    catch
        error('mstir library failed to load\n')
    end
end
end