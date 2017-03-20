function libload_stir
% load C++-to-C interface library
if ~libisloaded('mutilities')
    fprintf('loading mutilities library...\n')
    [notfound, warnings] = loadlibrary('mutilities');
end
% load STIR interface library
if ~libisloaded('mstir')
    fprintf('loading mstir library...\n')
    [notfound, warnings] = loadlibrary('mstir');
end
end