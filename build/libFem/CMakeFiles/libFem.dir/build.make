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
include libFem/CMakeFiles/libFem.dir/depend.make

# Include the progress variables for this target.
include libFem/CMakeFiles/libFem.dir/progress.make

# Include the compile flags for this target's objects.
include libFem/CMakeFiles/libFem.dir/flags.make

libFem/CMakeFiles/libFem.dir/src/libFEM.cpp.o: libFem/CMakeFiles/libFem.dir/flags.make
libFem/CMakeFiles/libFem.dir/src/libFEM.cpp.o: ../libFem/src/libFEM.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/liangqx/linearFEMGPU/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object libFem/CMakeFiles/libFem.dir/src/libFEM.cpp.o"
	cd /home/liangqx/linearFEMGPU/build/libFem && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/libFem.dir/src/libFEM.cpp.o -c /home/liangqx/linearFEMGPU/libFem/src/libFEM.cpp

libFem/CMakeFiles/libFem.dir/src/libFEM.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/libFem.dir/src/libFEM.cpp.i"
	cd /home/liangqx/linearFEMGPU/build/libFem && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/liangqx/linearFEMGPU/libFem/src/libFEM.cpp > CMakeFiles/libFem.dir/src/libFEM.cpp.i

libFem/CMakeFiles/libFem.dir/src/libFEM.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/libFem.dir/src/libFEM.cpp.s"
	cd /home/liangqx/linearFEMGPU/build/libFem && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/liangqx/linearFEMGPU/libFem/src/libFEM.cpp -o CMakeFiles/libFem.dir/src/libFEM.cpp.s

# Object files for target libFem
libFem_OBJECTS = \
"CMakeFiles/libFem.dir/src/libFEM.cpp.o"

# External object files for target libFem
libFem_EXTERNAL_OBJECTS =

libFem/liblibFem.a: libFem/CMakeFiles/libFem.dir/src/libFEM.cpp.o
libFem/liblibFem.a: libFem/CMakeFiles/libFem.dir/build.make
libFem/liblibFem.a: libFem/CMakeFiles/libFem.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/liangqx/linearFEMGPU/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library liblibFem.a"
	cd /home/liangqx/linearFEMGPU/build/libFem && $(CMAKE_COMMAND) -P CMakeFiles/libFem.dir/cmake_clean_target.cmake
	cd /home/liangqx/linearFEMGPU/build/libFem && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/libFem.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
libFem/CMakeFiles/libFem.dir/build: libFem/liblibFem.a

.PHONY : libFem/CMakeFiles/libFem.dir/build

libFem/CMakeFiles/libFem.dir/clean:
	cd /home/liangqx/linearFEMGPU/build/libFem && $(CMAKE_COMMAND) -P CMakeFiles/libFem.dir/cmake_clean.cmake
.PHONY : libFem/CMakeFiles/libFem.dir/clean

libFem/CMakeFiles/libFem.dir/depend:
	cd /home/liangqx/linearFEMGPU/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/liangqx/linearFEMGPU /home/liangqx/linearFEMGPU/libFem /home/liangqx/linearFEMGPU/build /home/liangqx/linearFEMGPU/build/libFem /home/liangqx/linearFEMGPU/build/libFem/CMakeFiles/libFem.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : libFem/CMakeFiles/libFem.dir/depend

