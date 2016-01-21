clear all

if libisloaded('mgadgetron')
    unloadlibrary('mgadgetron')
end

cgt_path = 'C:\Users\wps46139\Documents\GitHub\xGadgetron\cGadgetron';
cgt_lib = [cgt_path '\x64\Release\cGadgetron.lib'];

boost_ipath = getenv('BOOST');
boost_include = ['-I' boost_ipath];
boost_po_lib = 'C:\Boost\lib\libboost_program_options-vc120-mt-1_58.lib';
boost_system_lib = 'C:\Boost\lib\libboost_system-vc120-mt-1_58.lib';
boost_date_time_lib = 'C:\Boost\lib\libboost_date_time-vc120-mt-1_58.lib';
boost_regex_lib = 'C:\Boost\lib\libboost_regex-vc120-mt-1_58.lib';
boost_thread_lib = 'C:\Boost\lib\libboost_thread-vc120-mt-1_58.lib';
boost_chrono_lib = 'C:\Boost\lib\libboost_chrono-vc120-mt-1_58.lib';

ismrmrd_include = '-IC:\Users\wps46139\Documents\GitHub\ismrmrd\include';
ismrmrd_lib = 'C:\ISMRMRD\Release\ismrmrd.lib';

mex('-largeArrayDims',...
    boost_include, ismrmrd_include,...
    'mgadgetron.c', ...
    cgt_lib, ismrmrd_lib, ...
    boost_po_lib, boost_system_lib, boost_date_time_lib, ...
    boost_regex_lib, boost_thread_lib, boost_chrono_lib);
