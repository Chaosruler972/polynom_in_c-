# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.9

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
CMAKE_COMMAND = /home/chaosruler/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/173.3727.114/bin/cmake/bin/cmake

# The command to remove a file.
RM = /home/chaosruler/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/173.3727.114/bin/cmake/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/chaosruler/Desktop/polynom_in_c-

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/chaosruler/Desktop/polynom_in_c-/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/prepare_for_exam.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/prepare_for_exam.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/prepare_for_exam.dir/flags.make

CMakeFiles/prepare_for_exam.dir/main.cpp.o: CMakeFiles/prepare_for_exam.dir/flags.make
CMakeFiles/prepare_for_exam.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chaosruler/Desktop/polynom_in_c-/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/prepare_for_exam.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/prepare_for_exam.dir/main.cpp.o -c /home/chaosruler/Desktop/polynom_in_c-/main.cpp

CMakeFiles/prepare_for_exam.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/prepare_for_exam.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chaosruler/Desktop/polynom_in_c-/main.cpp > CMakeFiles/prepare_for_exam.dir/main.cpp.i

CMakeFiles/prepare_for_exam.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/prepare_for_exam.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chaosruler/Desktop/polynom_in_c-/main.cpp -o CMakeFiles/prepare_for_exam.dir/main.cpp.s

CMakeFiles/prepare_for_exam.dir/main.cpp.o.requires:

.PHONY : CMakeFiles/prepare_for_exam.dir/main.cpp.o.requires

CMakeFiles/prepare_for_exam.dir/main.cpp.o.provides: CMakeFiles/prepare_for_exam.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/prepare_for_exam.dir/build.make CMakeFiles/prepare_for_exam.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/prepare_for_exam.dir/main.cpp.o.provides

CMakeFiles/prepare_for_exam.dir/main.cpp.o.provides.build: CMakeFiles/prepare_for_exam.dir/main.cpp.o


CMakeFiles/prepare_for_exam.dir/polynom.cpp.o: CMakeFiles/prepare_for_exam.dir/flags.make
CMakeFiles/prepare_for_exam.dir/polynom.cpp.o: ../polynom.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chaosruler/Desktop/polynom_in_c-/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/prepare_for_exam.dir/polynom.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/prepare_for_exam.dir/polynom.cpp.o -c /home/chaosruler/Desktop/polynom_in_c-/polynom.cpp

CMakeFiles/prepare_for_exam.dir/polynom.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/prepare_for_exam.dir/polynom.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chaosruler/Desktop/polynom_in_c-/polynom.cpp > CMakeFiles/prepare_for_exam.dir/polynom.cpp.i

CMakeFiles/prepare_for_exam.dir/polynom.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/prepare_for_exam.dir/polynom.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chaosruler/Desktop/polynom_in_c-/polynom.cpp -o CMakeFiles/prepare_for_exam.dir/polynom.cpp.s

CMakeFiles/prepare_for_exam.dir/polynom.cpp.o.requires:

.PHONY : CMakeFiles/prepare_for_exam.dir/polynom.cpp.o.requires

CMakeFiles/prepare_for_exam.dir/polynom.cpp.o.provides: CMakeFiles/prepare_for_exam.dir/polynom.cpp.o.requires
	$(MAKE) -f CMakeFiles/prepare_for_exam.dir/build.make CMakeFiles/prepare_for_exam.dir/polynom.cpp.o.provides.build
.PHONY : CMakeFiles/prepare_for_exam.dir/polynom.cpp.o.provides

CMakeFiles/prepare_for_exam.dir/polynom.cpp.o.provides.build: CMakeFiles/prepare_for_exam.dir/polynom.cpp.o


# Object files for target prepare_for_exam
prepare_for_exam_OBJECTS = \
"CMakeFiles/prepare_for_exam.dir/main.cpp.o" \
"CMakeFiles/prepare_for_exam.dir/polynom.cpp.o"

# External object files for target prepare_for_exam
prepare_for_exam_EXTERNAL_OBJECTS =

prepare_for_exam: CMakeFiles/prepare_for_exam.dir/main.cpp.o
prepare_for_exam: CMakeFiles/prepare_for_exam.dir/polynom.cpp.o
prepare_for_exam: CMakeFiles/prepare_for_exam.dir/build.make
prepare_for_exam: CMakeFiles/prepare_for_exam.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/chaosruler/Desktop/polynom_in_c-/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable prepare_for_exam"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/prepare_for_exam.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/prepare_for_exam.dir/build: prepare_for_exam

.PHONY : CMakeFiles/prepare_for_exam.dir/build

CMakeFiles/prepare_for_exam.dir/requires: CMakeFiles/prepare_for_exam.dir/main.cpp.o.requires
CMakeFiles/prepare_for_exam.dir/requires: CMakeFiles/prepare_for_exam.dir/polynom.cpp.o.requires

.PHONY : CMakeFiles/prepare_for_exam.dir/requires

CMakeFiles/prepare_for_exam.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/prepare_for_exam.dir/cmake_clean.cmake
.PHONY : CMakeFiles/prepare_for_exam.dir/clean

CMakeFiles/prepare_for_exam.dir/depend:
	cd /home/chaosruler/Desktop/polynom_in_c-/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/chaosruler/Desktop/polynom_in_c- /home/chaosruler/Desktop/polynom_in_c- /home/chaosruler/Desktop/polynom_in_c-/cmake-build-debug /home/chaosruler/Desktop/polynom_in_c-/cmake-build-debug /home/chaosruler/Desktop/polynom_in_c-/cmake-build-debug/CMakeFiles/prepare_for_exam.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/prepare_for_exam.dir/depend

