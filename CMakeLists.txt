cmake_minimum_required(VERSION 3.16)

project(gravity LANGUAGES CXX)

# Various sets
set(CMAKE_CXX_STANDARD 17)

# Add Kokkos to the path
add_subdirectory(external/kokkos)

# Our sessions
add_subdirectory(session1)
add_subdirectory(session2)
add_subdirectory(session3)
add_subdirectory(session4)