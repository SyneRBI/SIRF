# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.7

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/sirfuser/devel/buildVM/sources/SIRF

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/sirfuser/devel/buildVM/sources/SIRF

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target install/strip
install/strip: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing the project stripped..."
	/usr/local/bin/cmake -DCMAKE_INSTALL_DO_STRIP=1 -P cmake_install.cmake
.PHONY : install/strip

# Special rule for the target install/strip
install/strip/fast: install/strip

.PHONY : install/strip/fast

# Special rule for the target install/local
install/local: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing only the local directory..."
	/usr/local/bin/cmake -DCMAKE_INSTALL_LOCAL_ONLY=1 -P cmake_install.cmake
.PHONY : install/local

# Special rule for the target install/local
install/local/fast: install/local

.PHONY : install/local/fast

# Special rule for the target test
test:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running tests..."
	/usr/local/bin/ctest --force-new-ctest-process $(ARGS)
.PHONY : test

# Special rule for the target test
test/fast: test

.PHONY : test/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/local/bin/ccmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/local/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target install
install: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Install the project..."
	/usr/local/bin/cmake -P cmake_install.cmake
.PHONY : install

# Special rule for the target install
install/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Install the project..."
	/usr/local/bin/cmake -P cmake_install.cmake
.PHONY : install/fast

# Special rule for the target list_install_components
list_install_components:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Available install components are: \"Unspecified\""
.PHONY : list_install_components

# Special rule for the target list_install_components
list_install_components/fast: list_install_components

.PHONY : list_install_components/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/sirfuser/devel/buildVM/sources/SIRF/CMakeFiles /home/sirfuser/devel/buildVM/sources/SIRF/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/sirfuser/devel/buildVM/sources/SIRF/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named iutilities

# Build rule for target.
iutilities: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 iutilities
.PHONY : iutilities

# fast build rule for target.
iutilities/fast:
	$(MAKE) -f src/iUtilities/CMakeFiles/iutilities.dir/build.make src/iUtilities/CMakeFiles/iutilities.dir/build
.PHONY : iutilities/fast

#=============================================================================
# Target rules for targets named _pyiutilities

# Build rule for target.
_pyiutilities: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 _pyiutilities
.PHONY : _pyiutilities

# fast build rule for target.
_pyiutilities/fast:
	$(MAKE) -f src/iUtilities/CMakeFiles/_pyiutilities.dir/build.make src/iUtilities/CMakeFiles/_pyiutilities.dir/build
.PHONY : _pyiutilities/fast

#=============================================================================
# Target rules for targets named cgadgetron

# Build rule for target.
cgadgetron: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 cgadgetron
.PHONY : cgadgetron

# fast build rule for target.
cgadgetron/fast:
	$(MAKE) -f src/xGadgetron/cGadgetron/CMakeFiles/cgadgetron.dir/build.make src/xGadgetron/cGadgetron/CMakeFiles/cgadgetron.dir/build
.PHONY : cgadgetron/fast

#=============================================================================
# Target rules for targets named _pygadgetron

# Build rule for target.
_pygadgetron: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 _pygadgetron
.PHONY : _pygadgetron

# fast build rule for target.
_pygadgetron/fast:
	$(MAKE) -f src/xGadgetron/pGadgetron/CMakeFiles/_pygadgetron.dir/build.make src/xGadgetron/pGadgetron/CMakeFiles/_pygadgetron.dir/build
.PHONY : _pygadgetron/fast

#=============================================================================
# Target rules for targets named cdynamicsimulation

# Build rule for target.
cdynamicsimulation: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 cdynamicsimulation
.PHONY : cdynamicsimulation

# fast build rule for target.
cdynamicsimulation/fast:
	$(MAKE) -f src/xDynamicSimulation/cDynamicSimulation/CMakeFiles/cdynamicsimulation.dir/build.make src/xDynamicSimulation/cDynamicSimulation/CMakeFiles/cdynamicsimulation.dir/build
.PHONY : cdynamicsimulation/fast

#=============================================================================
# Target rules for targets named TestJobs

# Build rule for target.
TestJobs: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 TestJobs
.PHONY : TestJobs

# fast build rule for target.
TestJobs/fast:
	$(MAKE) -f src/xDynamicSimulation/cDynamicSimulation/CMakeFiles/TestJobs.dir/build.make src/xDynamicSimulation/cDynamicSimulation/CMakeFiles/TestJobs.dir/build
.PHONY : TestJobs/fast

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... install/strip"
	@echo "... install/local"
	@echo "... test"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... install"
	@echo "... list_install_components"
	@echo "... iutilities"
	@echo "... _pyiutilities"
	@echo "... cgadgetron"
	@echo "... _pygadgetron"
	@echo "... cdynamicsimulation"
	@echo "... TestJobs"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

