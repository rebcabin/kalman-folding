# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.5.1/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.5.1/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/bbeckman/Documents/kalman-folding/C_CODE_PAPER_9/recurrenceForTheMean

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/bbeckman/Documents/kalman-folding/C_CODE_PAPER_9/recurrenceForTheMean

# Include any dependencies generated for this target.
include CMakeFiles/recurrenceForTheMean.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/recurrenceForTheMean.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/recurrenceForTheMean.dir/flags.make

CMakeFiles/recurrenceForTheMean.dir/main.c.o: CMakeFiles/recurrenceForTheMean.dir/flags.make
CMakeFiles/recurrenceForTheMean.dir/main.c.o: main.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/bbeckman/Documents/kalman-folding/C_CODE_PAPER_9/recurrenceForTheMean/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/recurrenceForTheMean.dir/main.c.o"
	/Library/Developer/CommandLineTools/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/recurrenceForTheMean.dir/main.c.o   -c /Users/bbeckman/Documents/kalman-folding/C_CODE_PAPER_9/recurrenceForTheMean/main.c

CMakeFiles/recurrenceForTheMean.dir/main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/recurrenceForTheMean.dir/main.c.i"
	/Library/Developer/CommandLineTools/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/bbeckman/Documents/kalman-folding/C_CODE_PAPER_9/recurrenceForTheMean/main.c > CMakeFiles/recurrenceForTheMean.dir/main.c.i

CMakeFiles/recurrenceForTheMean.dir/main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/recurrenceForTheMean.dir/main.c.s"
	/Library/Developer/CommandLineTools/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/bbeckman/Documents/kalman-folding/C_CODE_PAPER_9/recurrenceForTheMean/main.c -o CMakeFiles/recurrenceForTheMean.dir/main.c.s

CMakeFiles/recurrenceForTheMean.dir/main.c.o.requires:

.PHONY : CMakeFiles/recurrenceForTheMean.dir/main.c.o.requires

CMakeFiles/recurrenceForTheMean.dir/main.c.o.provides: CMakeFiles/recurrenceForTheMean.dir/main.c.o.requires
	$(MAKE) -f CMakeFiles/recurrenceForTheMean.dir/build.make CMakeFiles/recurrenceForTheMean.dir/main.c.o.provides.build
.PHONY : CMakeFiles/recurrenceForTheMean.dir/main.c.o.provides

CMakeFiles/recurrenceForTheMean.dir/main.c.o.provides.build: CMakeFiles/recurrenceForTheMean.dir/main.c.o


# Object files for target recurrenceForTheMean
recurrenceForTheMean_OBJECTS = \
"CMakeFiles/recurrenceForTheMean.dir/main.c.o"

# External object files for target recurrenceForTheMean
recurrenceForTheMean_EXTERNAL_OBJECTS =

recurrenceForTheMean: CMakeFiles/recurrenceForTheMean.dir/main.c.o
recurrenceForTheMean: CMakeFiles/recurrenceForTheMean.dir/build.make
recurrenceForTheMean: CMakeFiles/recurrenceForTheMean.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/bbeckman/Documents/kalman-folding/C_CODE_PAPER_9/recurrenceForTheMean/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable recurrenceForTheMean"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/recurrenceForTheMean.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/recurrenceForTheMean.dir/build: recurrenceForTheMean

.PHONY : CMakeFiles/recurrenceForTheMean.dir/build

CMakeFiles/recurrenceForTheMean.dir/requires: CMakeFiles/recurrenceForTheMean.dir/main.c.o.requires

.PHONY : CMakeFiles/recurrenceForTheMean.dir/requires

CMakeFiles/recurrenceForTheMean.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/recurrenceForTheMean.dir/cmake_clean.cmake
.PHONY : CMakeFiles/recurrenceForTheMean.dir/clean

CMakeFiles/recurrenceForTheMean.dir/depend:
	cd /Users/bbeckman/Documents/kalman-folding/C_CODE_PAPER_9/recurrenceForTheMean && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/bbeckman/Documents/kalman-folding/C_CODE_PAPER_9/recurrenceForTheMean /Users/bbeckman/Documents/kalman-folding/C_CODE_PAPER_9/recurrenceForTheMean /Users/bbeckman/Documents/kalman-folding/C_CODE_PAPER_9/recurrenceForTheMean /Users/bbeckman/Documents/kalman-folding/C_CODE_PAPER_9/recurrenceForTheMean /Users/bbeckman/Documents/kalman-folding/C_CODE_PAPER_9/recurrenceForTheMean/CMakeFiles/recurrenceForTheMean.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/recurrenceForTheMean.dir/depend
