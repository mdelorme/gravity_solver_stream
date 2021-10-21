add_executable(session6
  src/tree.cpp
  src/particles.cpp
  src/main.cpp)

target_link_libraries(session6 Kokkos::kokkos)