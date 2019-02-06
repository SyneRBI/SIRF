function libload_reg
% load C++-to-C interface library
if ~libisloaded('miutilities')
    fprintf('loading miutilities library...\n')
    loadlibrary('miutilities');
end
% load Reg interface library
if ~libisloaded('mreg')
    fprintf('loading mreg library...\n')
    loadlibrary('mreg');
end
end