cmake_minimum_required(VERSION 3.16)
project(ROOT_test
        VERSION 0.1
        DESCRIPTION "Event generator based on TFoam for DVCS on nuclei"
        HOMEPAGE_URL "https://gitlab.in2p3.fr/dupre/nopeg"
        LANGUAGES CXX Fortran
        )

# To force build in an other directory
if (CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR
            "In-source builds are not supported. "
            "Please read the BUILDING document before trying to build this project. "
            "You may need to delete 'CMakeCache.txt' and 'CMakeFiles/' first."
            )
endif ()

find_package(ROOT REQUIRED COMPONENTS RIO Net Hist Tree Matrix Physics MathCore Thread Foam GenVector)
include(${ROOT_USE_FILE})
set(SOURCE
        ${CMAKE_CURRENT_SOURCE_DIR}/source/Model131V2.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/source/Model531V2.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/source/TFoamMT.cpp
        )

set(HEADERS
        ${CMAKE_CURRENT_SOURCE_DIR}/include/Constant.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/draw.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/IO.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/Model.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/Model531V2.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/Parser.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/root2lund.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/TFoamMT.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/toml.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/Model131V2.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/BS_thread_pool_light.hpp
        )

set(FORTRAN_SOURCE
        ${CMAKE_CURRENT_SOURCE_DIR}/fortran/incoherent_deuteron_full_diff.f90
        )

find_package(fmt CONFIG REQUIRED)
add_executable(${PROJECT_NAME} source/main.cpp ${SOURCE} ${FORTRAN_SOURCE} ${HEADERS})

target_link_libraries(${PROJECT_NAME} PUBLIC ROOT::Core ROOT::Gpad ROOT::Physics ROOT::RIO ROOT::Tree Foam GenVector)
target_link_libraries(${PROJECT_NAME} PRIVATE fmt::fmt)
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fortran)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/source)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

message(STATUS "Configuration \n
System          ${CMAKE_SYSTEM}
Processor       ${PROCESSOR} (${CMAKE_SYSTEM_PROCESSOR})
Build type      ${CMAKE_BUILD_TYPE}
Install path    ${CMAKE_INSTALL_PREFIX}
Compiler        ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}
Compiler flags:
C               ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${_BUILD_TYPE_UPPER}}
C++             ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${_BUILD_TYPE_UPPER}}"
        )


