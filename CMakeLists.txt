cmake_minimum_required(VERSION 3.10)
project(ParticleLife)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Manually find system frameworks
set(SYSTEM_FRAMEWORK_PATH "/Library/Frameworks")

find_library(SFML_GRAPHICS_LIB NAMES sfml-graphics PATHS ${SYSTEM_FRAMEWORK_PATH} NO_DEFAULT_PATH)
find_library(SFML_WINDOW_LIB NAMES sfml-window PATHS ${SYSTEM_FRAMEWORK_PATH} NO_DEFAULT_PATH)
find_library(SFML_SYSTEM_LIB NAMES sfml-system PATHS ${SYSTEM_FRAMEWORK_PATH} NO_DEFAULT_PATH)
find_library(SFML_AUDIO_LIB NAMES sfml-audio PATHS ${SYSTEM_FRAMEWORK_PATH} NO_DEFAULT_PATH)
find_library(SFML_NETWORK_LIB NAMES sfml-network PATHS ${SYSTEM_FRAMEWORK_PATH} NO_DEFAULT_PATH)

if (NOT SFML_GRAPHICS_LIB OR NOT SFML_WINDOW_LIB OR NOT SFML_SYSTEM_LIB)
    message(FATAL_ERROR "Could not find SFML frameworks in ${SYSTEM_FRAMEWORK_PATH}")
endif()

message(STATUS "Found SFML Graphics: ${SFML_GRAPHICS_LIB}")

# Include directories
# Compilation worked previously, so standard paths might be enough.
# But we can add the framework path just in case.
set(CMAKE_FRAMEWORK_PATH ${SYSTEM_FRAMEWORK_PATH})

include_directories(particle)

# Source files
set(SOURCES
    particle/main.cpp
    particle/Universe.cpp
    particle/ResourcePath.cpp
)

# Executable
add_executable(ParticleLife ${SOURCES})

# Link SFML
target_link_libraries(ParticleLife
    ${SFML_GRAPHICS_LIB}
    ${SFML_WINDOW_LIB}
    ${SFML_SYSTEM_LIB}
    ${SFML_AUDIO_LIB}
    ${SFML_NETWORK_LIB}
)

# Add -F flag explicitly for the linker
target_link_options(ParticleLife PRIVATE "SHELL:-F${SYSTEM_FRAMEWORK_PATH}")
