myfilepath = mfilename('fullpath');
[mypath, myname, ext] = fileparts(myfilepath);
mUtilities.run_all_scripts(mypath)