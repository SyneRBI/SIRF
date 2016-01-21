function checkExecutionStatus(f, handle)
    status = calllib('mgadgetron', 'mExecutionStatus', handle);
    if status ~= 0
        errId = [f ':error'];
        msg = calllib('mgadgetron', 'mExecutionError', handle);
        file = calllib('mgadgetron', 'mExecutionErrorFile', handle);
        line = calllib('mgadgetron', 'mExecutionErrorLine', handle);
        error(errId, '%s at line %d of %s', msg, line, file)
    end
end
