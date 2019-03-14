myfilepath = mfilename('fullpath');
[mypath, myname, ext] = fileparts(myfilepath);
sirf.Utilities.run_all_scripts(mypath)