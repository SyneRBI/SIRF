function import_str = setup_MR(engine)
    if isempty(engine)
        engine = 'Gadgetron';
    end
    try
        eval(['libload_' lower(engine)])
    catch
        error('package %s not found\n', engine)
    end
    import_str = ['import m' engine '.*'];
end