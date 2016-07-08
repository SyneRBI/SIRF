function checkExecutionStatus(f, handle)
    status = calllib('mutilities', 'mExecutionStatus', handle);
    if status ~= 0
        errId = [f ':error'];
        msg = calllib('mutilities', 'mExecutionError', handle);
        file = calllib('mutilities', 'mExecutionErrorFile', handle);
        line = calllib('mutilities', 'mExecutionErrorLine', handle);
        error(errId, '%s at line %d of %s', msg, line, file)
    end
end
