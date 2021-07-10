add_executable(session5
  src/tree.cpp
  src/particles.cpp
  src/main.cpp)

target_link_libraries(session5 Kokkos::kokkos)