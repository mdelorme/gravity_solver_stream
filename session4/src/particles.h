#pragma once

#include <bits/stdc++.h>

#include <Kokkos_Core.hpp>

namespace gravity {

//!< Scratch pad memory array
using ScratchPadView 
  = Kokkos::View<double**, 
                 Kokkos::DefaultExecutionSpace::scratch_memory_space, 
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

// Particle table : double[N][10]
enum Fields : uint8_t {
  IX = 0,
  IY = 1,
  IZ = 2,
  IVX = 3,
  IVY = 4,
  IVZ = 5,
  IAX = 6,
  IAY = 7,
  IAZ = 8,
  IM  = 9,

  FIELD_COUNT=10
};

using DataArray     = Kokkos::View<double**, Kokkos::LayoutRight>;
using DataArrayHost = DataArray::HostMirror;
using PartHost = std::array<double, 10>;

struct Particles {
  double G;
  DataArray     data;   //!< Données en mémoire sur le GPU
  DataArrayHost data_h; //!< Données en mémoire sur le CPU
  uint N;

  Particles(std::vector<PartHost> &input);

  void Send2GPU();
  void Send2CPU();

  void ResetAccelerations();
  void ComputeAccelerations();
  void Update(double dt);
};


}