# CMake generated Testfile for 
# Source directory: /home/sirfuser/devel/buildVM/sources/SIRF/src/xGadgetron/pGadgetron/tests
# Build directory: /home/sirfuser/devel/buildVM/sources/SIRF/src/xGadgetron/pGadgetron/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(MR_TESTS_PYTHON "/usr/bin/python" "-m" "nose" "--with-coverage" "--cover-package=pGadgetron" "src/xGadgetron/pGadgetron/")
set_tests_properties(MR_TESTS_PYTHON PROPERTIES  WORKING_DIRECTORY "/home/sirfuser/devel/buildVM/sources/SIRF")
