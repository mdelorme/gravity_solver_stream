#pragma once

#include <Kokkos_Core.hpp>
#include <bits/stdc++.h>

namespace gravity {

constexpr double theta_crit = 1.0;

// Particle table : double[N][10]
enum Fields : uint8_t {
  IX  = 0,
  IY  = 1,
  IZ  = 2,
  IM  = 3,
  IVX = 4,
  IVY = 5,
  IVZ = 6,
  IAX = 7,
  IAY = 8,
  IAZ = 9,

  FIELD_COUNT=10
};

using DataArray = Kokkos::View<double**>;
using PartHost  = std::array<double, FIELD_COUNT>;
using TreePart  = std::array<double, 4>;
using Vec3      = std::array<double, 3>;

}