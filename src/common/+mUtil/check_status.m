function check_status(f, handle)
    lib = 'mutilities';
    status = calllib(lib, 'mExecutionStatus', handle);
    if status ~= 0
        errId = [f ':error'];
        msg = calllib(lib, 'mExecutionError', handle);
        file = calllib(lib, 'mExecutionErrorFile', handle);
        line = calllib(lib', 'mExecutionErrorLine', handle);
        error(errId, '??? %s exception thrown at line %d of %s', ...
            msg, line, file)
    end
end
