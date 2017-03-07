function select_gadgetron
% Selects Gadgetron as MR Engine

filepath = mfilename('fullpath');
l = length(filepath) - length('select_gadgetron');
path = filepath(1:l);
copyfile([path '/+mGadgetron'], [path '/+MR'], 'f')

if ~libisloaded('mutilities')
    fprintf('loading mutilities library...\n')
    [notfound, warnings] = loadlibrary('mutilities');
end
if ~libisloaded('mgadgetron')
    fprintf('loading mgadgetron library...\n')
    [notfound, warnings] = loadlibrary('mgadgetron');
end

%libfunctions('mutilities')
%libfunctions('mgadgetron')
end