myfilepath = mfilename('fullpath');
[mypath, myname, ext] = fileparts(myfilepath);
run_all_scripts(mypath)