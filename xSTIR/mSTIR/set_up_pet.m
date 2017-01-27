try
    eval(['select_' PET_ENGINE])
catch
    if exist('PET_ENGINE', 'var')
        error('package %s not found\n', PET_ENGINE)
    else
        select_stir
    end
end
