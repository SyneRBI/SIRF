function ccp_libload
% CCP_LIBLOAD Load mutilities and mgadgetron libraries for SIRF/CCP code
%
% Usage:
%  ccp_libload
%
% Extracted by David Atkinson from code by Evgueni Ovtchinnikov
% 

% load the SIRF mutilities library
if ~libisloaded('mutilities')
    fprintf('loading mutilities library...\n')
    [notfound, warnings] = loadlibrary('mutilities');
end

 
% Load the SIRF library for gadgetron
if ~libisloaded('mgadgetron')
    fprintf('loading mgadgetron library...\n')
    [notfound, warnings] = loadlibrary('mgadgetron');
end