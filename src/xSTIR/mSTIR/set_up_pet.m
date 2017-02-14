% Sets up PET environment.
% For quick access to data files assign path to SIRF folder to SIRF_PATH.

if ~exist('SIRF_PET_DATA_PATH', 'var')
    if exist('SIRF_PATH', 'var')
        SIRF_PET_DATA_PATH = [SIRF_PATH '/data/examples/PET'];
    else
        SIRF_PET_DATA_PATH = './';
    end
end

try
    eval(['select_' PET_ENGINE])
catch
    if exist('PET_ENGINE', 'var')
        error('package %s not found\n', PET_ENGINE)
    else
        select_stir
    end
end
