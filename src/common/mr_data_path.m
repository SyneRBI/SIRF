function p = mr_data_path
% Tries to find path to MR raw data.
% The user may like to set a Matlab variable SIRF_MR_DATA_PATH
% to the path to their raw MR data.
% If it is not set, the path to SIRF subfolder /data/examples/MR
% will be used.

if exist('SIRF_MR_DATA_PATH', 'var')
    p = SIRF_MR_DATA_PATH;
else
    SIRF_PATH = getenv('SIRF_PATH');
    if ~isempty(SIRF_PATH)
        p = [SIRF_PATH '/data/examples/MR'];
    else
        p = './';
    end
end
end