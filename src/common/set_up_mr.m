function set_up_mr
% Sets up MR Engine

try
    eval(['select_' MR_ENGINE])
catch
    if exist('MR_ENGINE', 'var')
        error('package %s not found\n', MR_ENGINE)
    else
        select_gadgetron
    end
end
end