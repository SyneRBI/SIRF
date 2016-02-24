clear all

if libisloaded('mgadgetron')
    unloadlibrary('mgadgetron')
end
if libisloaded('mutilities')
    unloadlibrary('mutilities')
end

cgt_path = getenv('CGADGETRON');
cgt_lib = getenv('CGADGETRON_LIBRARY');

util_path = getenv('IUTILITIES');
util_lib = getenv('IUTILITIES_LIBRARY');
util_include = ['-I' util_path];
tw = [util_path '/text_writer.cpp'];

boost_ipath = getenv('BOOST');
boost_lpath = getenv('BOOST_LIB');
boost_suffix = getenv('BOOST_SUFFIX');
boost_include = ['-I' boost_ipath];
boost_po_lib = [boost_lpath '/libboost_program_options' boost_suffix];
boost_system_lib = [boost_lpath '/libboost_system' boost_suffix];
boost_date_time_lib = [boost_lpath '/libboost_date_time' boost_suffix];
boost_regex_lib = [boost_lpath '/libboost_regex' boost_suffix];
boost_thread_lib = [boost_lpath '/libboost_thread' boost_suffix];
boost_chrono_lib = [boost_lpath '/libboost_chrono' boost_suffix];

ismrmrd_include = ['-I' getenv('ISMRMRD_INCLUDE')];
ismrmrd_lib = [getenv('ISMRMRD_LIB') '/' getenv('ISMRMRD_LIBRARY')];

if strcmp(getenv('GCC'), 'gcc')
    CCFLAG = '-DGCC';
else
    CCFLAG = '-DNOGCC';
end

mex('-largeArrayDims', CCFLAG, ...
    boost_include, util_include, ...
    'mutilities.c', 'printer.cpp', ...
    util_lib) 

mex('-largeArrayDims', CCFLAG, ...
    boost_include, ismrmrd_include, util_include, ...
    'mgadgetron.c', tw, ...
    cgt_lib, util_lib, ismrmrd_lib, ...
    boost_po_lib, boost_system_lib, boost_date_time_lib, ...
    boost_regex_lib, boost_thread_lib, boost_chrono_lib);
