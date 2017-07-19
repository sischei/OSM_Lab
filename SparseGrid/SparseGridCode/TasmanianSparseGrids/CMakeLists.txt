cmake_minimum_required(VERSION 2.8)
project(TasmanianSparseGrids)

#find_package(CUDA QUIET REQUIRED)

# define source files and targets (libraries and executables)
add_executable(tasgrid tasgrid_main.cpp
                       TasmanianSparseGrid.hpp
                       tasgridExternalTests.hpp
                       tasgridExternalTests.cpp
                       tasgridTestFunctions.hpp
                       tasgridTestFunctions.cpp
                       tasgridWrapper.hpp
                       tasgridWrapper.cpp)

add_library(libtas_object OBJECT TasmanianSparseGrid.hpp
                                 TasmanianSparseGrid.cpp
                                 tsgCacheLagrange.hpp
                                 tsgCacheLagrange.cpp
                                 tsgCoreOneDimensional.hpp
                                 tsgCoreOneDimensional.cpp
                                 tsgEnumerates.hpp
                                 tsgGridCore.hpp
                                 tsgGridCore.cpp
                                 tsgGridGlobal.hpp
                                 tsgGridGlobal.cpp
                                 tsgGridWavelet.hpp
                                 tsgGridWavelet.cpp
                                 tsgHardCodedTabulatedRules.hpp
                                 tsgHardCodedTabulatedRules.cpp
                                 tsgGridLocalPolynomial.hpp
                                 tsgGridLocalPolynomial.cpp
                                 tsgGridSequence.hpp
                                 tsgGridSequence.cpp
                                 tsgIndexManipulator.hpp
                                 tsgIndexManipulator.cpp
                                 tsgIndexSets.hpp
                                 tsgIndexSets.cpp
                                 tsgLinearSolvers.hpp
                                 tsgLinearSolvers.cpp
                                 tsgOneDimensionalWrapper.hpp
                                 tsgOneDimensionalWrapper.cpp
                                 tsgRuleLocalPolynomial.hpp
                                 tsgRuleLocalPolynomial.cpp
                                 tsgRuleWavelet.hpp
                                 tsgRuleWavelet.cpp
                                 tsgSequenceOptimizer.hpp
                                 tsgSequenceOptimizer.cpp )

set_property(TARGET libtas_object PROPERTY POSITION_INDEPENDENT_CODE 1)

add_library(libtas_shared SHARED $<TARGET_OBJECTS:libtas_object>)
add_library(libtas_static STATIC $<TARGET_OBJECTS:libtas_object>)

target_link_libraries(tasgrid libtas_static)

# I like to keep the "lib" in front on a Windows machine
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    set(TSGBaseLibName "libtasmaniansparsegrid")
else()
    set(TSGBaseLibName "tasmaniansparsegrid")
endif()
set_target_properties(libtas_shared PROPERTIES OUTPUT_NAME "${TSGBaseLibName}")
set_target_properties(libtas_static PROPERTIES OUTPUT_NAME "${TSGBaseLibName}")

# options
option(ENABLE_OPENMP "ENABLE_OPENMP" ON)
option(ENABLE_MATLAB "ENABLE_MATLAB" ON)
option(ENABLE_PYTHON "ENABLE_PYTHON" ON)

if ( ENABLE_OPENMP )
    find_package(OpenMP)
    if (OPENMP_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
    endif()
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
endif()

# compile options
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ffast-math -mtune=native -Wall -Wextra -Wshadow -Wno-unused-parameter -pedantic")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mtune=native -diag-disable 11074 -diag-disable 11076 -Wall -Wextra -Wshadow -Wno-unused-parameter -pedantic")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ffast-math -mtune=native -Wall -Wextra -Wshadow -Wno-unused-parameter -pedantic")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Ox /EHsc -D_SCL_SECURE_NO_WARNINGS")
endif()

# paths, filenames, variables to add to .in. files
set(TSGLibPath "${CMAKE_INSTALL_PREFIX}/lib")
set(TSGLibShared "${TSGLibPath}/${CMAKE_SHARED_LIBRARY_PREFIX}${TSGBaseLibName}${CMAKE_SHARED_LIBRARY_SUFFIX}")
set(TSGLibStatic "${TSGLibPath}/${CMAKE_STATIC_LIBRARY_PREFIX}${TSGBaseLibName}${CMAKE_STATIC_LIBRARY_SUFFIX}")
set(TSGExecutable "${CMAKE_INSTALL_PREFIX}/bin/tasgrid")
set(TSGIncludePath "${CMAKE_INSTALL_PREFIX}/include")
set(TSGPythonPath "${CMAKE_INSTALL_PREFIX}/python")
set(TSGMatlabPath "${CMAKE_INSTALL_PREFIX}/matlab")
set(TSGExamplePath "${CMAKE_INSTALL_PREFIX}/example")
set(TSGLibOPTC "${CMAKE_CXX_FLAGS}")
set(TSGExampleImport "import sys\nsys.path.append(\"${TSGPythonPath}\")\n")
set(OMPBasicTest "")

if ( NOT DEFINED MATLAB_WORK_FOLDER )
    set(TSGMatlabWork "${CMAKE_INSTALL_PREFIX}/tsgMatlabWork/")
else()
    set(TSGMatlabWork "${MATLAB_WORK_FOLDER}")
endif()

# Gauss-Pattersion .table file and Python interface
configure_file( "${PROJECT_SOURCE_DIR}/GaussPattersonRule.table"  "${CMAKE_BINARY_DIR}/GaussPattersonRule.table" COPYONLY )
if ( ENABLE_PYTHON )
    configure_file( "${PROJECT_SOURCE_DIR}/InterfacePython/TasmanianSG.in.py"  "${CMAKE_BINARY_DIR}/python/TasmanianSG.py" )
    configure_file( "${PROJECT_SOURCE_DIR}/InterfacePython/example.in.py"  "${CMAKE_BINARY_DIR}/example/example.py" )
endif()
if ( ENABLE_MATLAB )
    configure_file( "${PROJECT_SOURCE_DIR}/InterfaceMATLAB/tsgGetPaths.in.m"  "${CMAKE_BINARY_DIR}/matlab/tsgGetPaths.m" )
endif()
configure_file( "${PROJECT_SOURCE_DIR}/TasmanianSparseGrid.h"  "${CMAKE_BINARY_DIR}/include/TasmanianSparseGrid.h" COPYONLY )
configure_file( "${PROJECT_SOURCE_DIR}/Example/Makefile.in"  "${CMAKE_BINARY_DIR}/example/Makefile" )
configure_file( "${PROJECT_SOURCE_DIR}/Example/example.cpp"  "${CMAKE_BINARY_DIR}/example/example.cpp" COPYONLY )


# for support of simple GNU-Make
#set(GENERATE_BACKWARDS "1")
if ( DEFINED GENERATE_BACKWARDS )
    set(TSGLibShared "../libtasmaniansparsegrid.so")
    configure_file( "${PROJECT_SOURCE_DIR}/InterfacePython/TasmanianSG.in.py"  "${PROJECT_SOURCE_DIR}/InterfacePython/TasmanianSG.py" )
    set(TSGExampleImport "#")
    configure_file( "${PROJECT_SOURCE_DIR}/InterfacePython/example.in.py"  "${PROJECT_SOURCE_DIR}/InterfacePython/example.py" )
    set(TSGExecutable "ENTER THE PATH TO tasgrid EXECUTABLE")
    set(TSGMatlabWork "ENTER THE PATH TO MATLAB WORK FOLDER")
    configure_file( "${PROJECT_SOURCE_DIR}/InterfaceMATLAB/tsgGetPaths.in.m"  "${PROJECT_SOURCE_DIR}/InterfaceMATLAB/tsgGetPaths.m" )
    set(OMPBasicTest "UNAME = \$(shell uname)\nifeq (\$(UNAME), Linux)\nCOMPILE_OPTIONS = -fopenmp\nelse\nCOMPILE_OPTIONS = -Wno-unknown-pragmas\nendif\n")
    set(CMAKE_CXX_COMPILER "g++")
    set(TSGLibOPTC "-O3 -ffast-math -mtune=native $(COMPILE_OPTIONS) -fPIC -Wall -Wextra -Wshadow -Wno-unused-parameter -pedantic")
    set(TSGLibStatic "../libtasmaniansparsegrid.a")
    set(TSGIncludePath "../")
    configure_file( "${PROJECT_SOURCE_DIR}/Example/Makefile.in"  "${PROJECT_SOURCE_DIR}/Example/Makefile" )
endif()

## for the purpose of development-debugging
set(TSGLibShared "./libtasmaniansparsegrid.so")
configure_file( "${PROJECT_SOURCE_DIR}/InterfacePython/TasmanianSG.in.py"  "${CMAKE_BINARY_DIR}/TasmanianSG.py" )
set(TSGLibShared "${TSGLibPath}/${CMAKE_SHARED_LIBRARY_PREFIX}${TSGBaseLibName}${CMAKE_SHARED_LIBRARY_SUFFIX}")
set(TSGExampleImport "#")
configure_file( "${PROJECT_SOURCE_DIR}/InterfacePython/example.in.py"  "${CMAKE_BINARY_DIR}/example.py" )
configure_file( "${PROJECT_SOURCE_DIR}/InterfacePython/sandbox.py"  "${CMAKE_BINARY_DIR}/sandbox.py" COPYONLY)
set(TSGExampleImport "import sys\nsys.path.append(\"${TSGPythonPath}\")\n")
configure_file( "${PROJECT_SOURCE_DIR}/InterfacePython/tsgTests.py"  "${CMAKE_BINARY_DIR}/tsgTests.py" COPYONLY)

enable_testing()
add_test(Tasgrid tasgrid -test)

# install commands
install(TARGETS tasgrid libtas_shared libtas_static
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/
        DESTINATION include
        FILES_MATCHING PATTERN "*.hpp"
        PATTERN "InterfacePython" EXCLUDE
        PATTERN "InterfaceMATLAB" EXCLUDE
        PATTERN "HelperScripts" EXCLUDE
        PATTERN "Example" EXCLUDE)
install(FILES ${CMAKE_BINARY_DIR}/GaussPattersonRule.table
        DESTINATION example
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)

# C++ example
install(FILES ${CMAKE_BINARY_DIR}/example/example.cpp
        DESTINATION example
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
install(FILES ${CMAKE_BINARY_DIR}/example/Makefile
        DESTINATION example
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)

# C Header
install(FILES ${CMAKE_SOURCE_DIR}/TasmanianSparseGrid.h
        DESTINATION include
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)

if ( ENABLE_PYTHON )
install(FILES ${CMAKE_BINARY_DIR}/python/TasmanianSG.py
        DESTINATION python
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
install(FILES ${CMAKE_BINARY_DIR}/example/example.py
        DESTINATION example
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
endif()

if ( ENABLE_MATLAB )
install(DIRECTORY ${CMAKE_SOURCE_DIR}/InterfaceMATLAB/
        DESTINATION matlab
        FILES_MATCHING PATTERN "*.m"
        PATTERN "tsgGetPaths.*" EXCLUDE)
install(FILES ${CMAKE_BINARY_DIR}/matlab/tsgGetPaths.m
        DESTINATION matlab
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
endif()

# final message
message("----------------------------------------------------------")
message("Tasmanian Sparse Grids setup complete\n")
message(" executable:     ${TSGExecutable}")
message(" shared library: ${TSGLibShared}")
message(" static library: ${TSGLibStatic}")
message(" C++ example:    ${CMAKE_BINARY_DIR}/example/example.cpp")
message("                 see \"Makefile\" command in ${CMAKE_BINARY_DIR}/example")
message("")
message(" bash export commands (add to your environment, i.e., ~/.bashrc):")
message("")
message("      export PATH=$PATH:${CMAKE_INSTALL_PREFIX}/bin/")
message("      export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CMAKE_INSTALL_PREFIX}/lib/")
message("      export C_INCLUDE_PATH=$C_INCLUDE_PATH:${CMAKE_INSTALL_PREFIX}/include/")
message("      export CPP_INCLUDE_PATH=$CPP_INCLUDE_PATH:${CMAKE_INSTALL_PREFIX}/include/")
message("")
message("")

if ( ENABLE_PYTHON )
    message(" python module:  ${CMAKE_INSTALL_PREFIX}/python/TasmanianSG.py")
    message(" python example: ${CMAKE_INSTALL_PREFIX}/example/example.py")
    message("                 to call the python module")
    message("")
    message("                 import sys")
    message("                 sys.path.append(\"${CMAKE_INSTALL_PREFIX}/python\")")
    message("                 import TasmanianSG")
    message("")
    message("")
endif()

if ( ENABLE_MATLAB )
    message(" matlab flles in: ${CMAKE_INSTALL_PREFIX}/matlab/")
    message("                  NOTE: add ${CMAKE_INSTALL_PREFIX}/matlab/ to the MATLAB path")
    message("")
    message("                        addpath('${CMAKE_INSTALL_PREFIX}/matlab/')")
    message("")
    message("")
endif()
message("----------------------------------------------------------")

