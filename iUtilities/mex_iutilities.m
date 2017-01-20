clear all

if libisloaded('mutilities')
    unloadlibrary('mutilities')
end

boost_ipath = getenv('BOOST');
boost_lpath = getenv('BOOST_LIB');
boost_suffix = getenv('BOOST_SUFFIX');
boost_include = ['-I' boost_ipath];

util_path = getenv('IUTILITIES');
util_lib = getenv('IUTILITIES_LIBRARY');
util_include = ['-I' util_path];

mex('-largeArrayDims', ...
    boost_include, util_include, ...
    'mutilities.c', ...
    util_lib) 
