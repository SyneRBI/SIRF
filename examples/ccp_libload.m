function ccp_libload
% CCP_LIBLOAD Load mutilities and mgadgetron libraries for CCP code
%
% Usage:
%  ccp_libload
%
% Extracted by David Atkinson from code by Evgueni Ovtchinnikov
% 

% mutilities is now in folder iUtilities
if ~libisloaded('mutilities')
    fprintf('loading mutilities library...\n')
    [notfound, warnings] = loadlibrary('mutilities');
end

if ~libisloaded('mgadgetron')
    fprintf('loading mgadgetron library...\n')
    [notfound, warnings] = loadlibrary('mgadgetron');
end