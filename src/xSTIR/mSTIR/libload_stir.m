function libload_stir
% load C++-to-C interface library
if ~libisloaded('miutilities')
    fprintf('loading miutilities library...\n')
    [notfound, warnings] = loadlibrary('miutilities');
end
% load STIR interface library
if ~libisloaded('mstir')
    fprintf('loading mstir library...\n')
    [notfound, warnings] = loadlibrary('mstir');
end
end