function libload_gadgetron
% Loads Gadgetron interface libraries

% if ~libisloaded('mutilities')
%     fprintf('loading mutilities library...\n')
%     [notfound, warnings] = loadlibrary('mutilities');
% end
if ~libisloaded('miutilities')
    fprintf('loading miutilities library...\n')
    loadlibrary('miutilities');
end
if ~libisloaded('mgadgetron')
    fprintf('loading mgadgetron library...\n')
    loadlibrary('mgadgetron');
end
end