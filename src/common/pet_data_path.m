function p = pet_data_path
% For quick access to data files assign path to SIRF folder to SIRF_PATH.

if exist('SIRF_PET_DATA_PATH', 'var')
    p = SIRF_PET_DATA_PATH;
else
    if exist('SIRF_PATH', 'var')
        p = [SIRF_PATH '/data/examples/PET'];
    else
        p = './';
    end
end
end