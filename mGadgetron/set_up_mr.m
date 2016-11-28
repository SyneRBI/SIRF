try
    eval(['select_' MR_ENGINE])
catch
    if exist('MR_ENGINE', 'var')
        error('package %s not found\n', MR_ENGINE)
    else
        import gadgetron.*
    end
end
