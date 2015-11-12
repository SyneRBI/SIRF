function checkExecutionStatus(f, handle)
    status = calllib('mstir', 'mExecutionStatus', handle);
    if status ~= 0
        errId = [f ':error'];
        msg = calllib('mstir', 'mExecutionError', handle);
        file = calllib('mstir', 'mExecutionErrorFile', handle);
        line = calllib('mstir', 'mExecutionErrorLine', handle);
        error(errId, '%s at line %d of %s\n', msg, line, file);
    end
end
