myfilepath = mfilename('fullpath');
[mypath, myname, ext] = fileparts(myfilepath);
filelist = dir(fullfile(mypath, '*.m'));
nf = length(filelist);
failed = 0;
ntests = 0;
for k = 1 : nf
    filename = filelist(k).name;
    if ~strcmpi(filename, [myname, '.m'])
        [path, name, ext] = fileparts(filename);
        fprintf('Running %s\n', filename)
        [f, n] = eval(name);
        failed = failed + f;
        ntests = ntests + n;
        fprintf('-----------------\n')
    end
end
if failed == 0
    fprintf('all %d tests passed\n', ntests);
else
    fprintf('%d out of %d tests failed\n', failed, ntests);
end
