function libload_sirf
% load C++-to-C interface library
if ~libisloaded('miutilities')
    fprintf('loading miutilities library...\n')
    loadlibrary('miutilities');
end
% load STIR interface library
if ~libisloaded('msirf')
    fprintf('loading msirf library...\n')
    loadlibrary('msirf');
end
end