# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/liangqx/linearFEMGPU

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/liangqx/linearFEMGPU/build

# Include any dependencies generated for this target.
include preProcessor/CMakeFiles/preprocessor.dir/depend.make

# Include the progress variables for this target.
include preProcessor/CMakeFiles/preprocessor.dir/progress.make

# Include the compile flags for this target's objects.
include preProcessor/CMakeFiles/preprocessor.dir/flags.make

preProcessor/CMakeFiles/preprocessor.dir/src/preprocessor.cpp.o: preProcessor/CMakeFiles/preprocessor.dir/flags.make
preProcessor/CMakeFiles/preprocessor.dir/src/preprocessor.cpp.o: ../preProcessor/src/preprocessor.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/liangqx/linearFEMGPU/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object preProcessor/CMakeFiles/preprocessor.dir/src/preprocessor.cpp.o"
	cd /home/liangqx/linearFEMGPU/build/preProcessor && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/preprocessor.dir/src/preprocessor.cpp.o -c /home/liangqx/linearFEMGPU/preProcessor/src/preprocessor.cpp

preProcessor/CMakeFiles/preprocessor.dir/src/preprocessor.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/preprocessor.dir/src/preprocessor.cpp.i"
	cd /home/liangqx/linearFEMGPU/build/preProcessor && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/liangqx/linearFEMGPU/preProcessor/src/preprocessor.cpp > CMakeFiles/preprocessor.dir/src/preprocessor.cpp.i

preProcessor/CMakeFiles/preprocessor.dir/src/preprocessor.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/preprocessor.dir/src/preprocessor.cpp.s"
	cd /home/liangqx/linearFEMGPU/build/preProcessor && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/liangqx/linearFEMGPU/preProcessor/src/preprocessor.cpp -o CMakeFiles/preprocessor.dir/src/preprocessor.cpp.s

# Object files for target preprocessor
preprocessor_OBJECTS = \
"CMakeFiles/preprocessor.dir/src/preprocessor.cpp.o"

# External object files for target preprocessor
preprocessor_EXTERNAL_OBJECTS =

preProcessor/libpreprocessor.a: preProcessor/CMakeFiles/preprocessor.dir/src/preprocessor.cpp.o
preProcessor/libpreprocessor.a: preProcessor/CMakeFiles/preprocessor.dir/build.make
preProcessor/libpreprocessor.a: preProcessor/CMakeFiles/preprocessor.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/liangqx/linearFEMGPU/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libpreprocessor.a"
	cd /home/liangqx/linearFEMGPU/build/preProcessor && $(CMAKE_COMMAND) -P CMakeFiles/preprocessor.dir/cmake_clean_target.cmake
	cd /home/liangqx/linearFEMGPU/build/preProcessor && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/preprocessor.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
preProcessor/CMakeFiles/preprocessor.dir/build: preProcessor/libpreprocessor.a

.PHONY : preProcessor/CMakeFiles/preprocessor.dir/build

preProcessor/CMakeFiles/preprocessor.dir/clean:
	cd /home/liangqx/linearFEMGPU/build/preProcessor && $(CMAKE_COMMAND) -P CMakeFiles/preprocessor.dir/cmake_clean.cmake
.PHONY : preProcessor/CMakeFiles/preprocessor.dir/clean

preProcessor/CMakeFiles/preprocessor.dir/depend:
	cd /home/liangqx/linearFEMGPU/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/liangqx/linearFEMGPU /home/liangqx/linearFEMGPU/preProcessor /home/liangqx/linearFEMGPU/build /home/liangqx/linearFEMGPU/build/preProcessor /home/liangqx/linearFEMGPU/build/preProcessor/CMakeFiles/preprocessor.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : preProcessor/CMakeFiles/preprocessor.dir/depend

