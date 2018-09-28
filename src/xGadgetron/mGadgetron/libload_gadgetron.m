function libload_gadgetron
% Loads Gadgetron interface libraries

% if ~libisloaded('mutilities')
%     fprintf('loading mutilities library...\n')
%     [notfound, warnings] = loadlibrary('mutilities');
% end
if ~libisloaded('miutilities')
    fprintf('loading miutilities library...\n')
    try
        [notfound, warnings] = loadlibrary('miutilities');
    catch
        error('miutilities library failed to load\n')
    end
end
if ~libisloaded('mgadgetron')
    fprintf('loading mgadgetron library...\n')
    try
        [notfound, warnings] = loadlibrary('mgadgetron');
    catch
        error('mgadgetron library failed to load\n')
    end
end
end