#pragma once

#include <bits/stdc++.h>
#include "tree.h"
#include "params.h"

namespace gravity {

struct Particles {
  double G;
  DataArray data;   //!< Données en mémoire sur le CPU
  uint N;

  Particles(std::vector<PartHost> &input);

  void ResetAccelerations();
  void ComputeAccelerations(const Cell &tree);
  void Update(double dt);
};


}