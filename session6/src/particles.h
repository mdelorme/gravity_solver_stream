#pragma once

#include <bits/stdc++.h>
#include "tree.h"
#include "params.h"

namespace gravity {

struct Particles {
  double G;
  DataArray data;   //!< Données en mémoire sur le CPU
  Params params;    //!< Paramètres du run
  uint N;           //!< Nombre de particules

  Particles(std::vector<PartHost> &input, const Params& params);

  void ResetAccelerations();
  void ComputeAccelerations(const Cell &tree);
  void Update(double dt);
};


}