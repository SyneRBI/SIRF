function run_all_scripts(path)
filelist = dir(fullfile(path, '*.m'));
nf = length(filelist);
for k = 1 : nf
    filename = filelist(k).name;
    if ~strcmpi(filename, 'run_all.m')
        fprintf('Running %s\n', filename)
        try
            run(fullfile(path, filename));
            fprintf('ok\n')
        catch
            fprintf('failed\n')
        end
    end
end