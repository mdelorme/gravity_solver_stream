#include "particles.h"

namespace gravity {

/**
 * Foncteur 1: Mise à zero des accélérations
 **/
struct ResetAccelerationsFunctor {
  ResetAccelerationsFunctor(DataArray d) : data(d) {};
  
  KOKKOS_INLINE_FUNCTION
  void operator()(const uint i) const {
    data(i, IAX) = 0.0;
    data(i, IAY) = 0.0;
    data(i, IAZ) = 0.0;    
  }

  static void apply(DataArray d, uint N) {
    ResetAccelerationsFunctor functor(d);
    Kokkos::RangePolicy policy(0, N);
    Kokkos::parallel_for(policy, functor);
  }

  // Public members
  DataArray data;
};

/**
 * Foncteur 2: Calcul des accélerations
 **/
struct ComputeAccelerationsFunctor {
  ComputeAccelerationsFunctor(DataArray d, const Cell &tree, 
                              double G, double theta_crit) 
    : data(d), tree(tree), G(G), theta_crit(theta_crit) {};

  KOKKOS_INLINE_FUNCTION
  void operator()(const uint i) const {
    Vec3 xi;
    xi[IX] = data(i, IX);
    xi[IY] = data(i, IY);
    xi[IZ] = data(i, IZ);

    Vec3 ai = tree.CalculateAcceleration(xi, G, theta_crit);
    data(i, IAX) = ai[IX];
    data(i, IAY) = ai[IY];
    data(i, IAZ) = ai[IZ];
  }

  static void apply(DataArray d, const Cell &tree, uint N, 
                    double G, double theta_crit) {
    ComputeAccelerationsFunctor functor(d, tree, G, theta_crit);
    Kokkos::RangePolicy policy(0, N);
    Kokkos::parallel_for(policy, functor);
  }

  // Public members
  DataArray   data;  
  const Cell &tree;
  double G;
  double theta_crit;
};

/**
 * Foncteur 3: Mise à jour
 **/
struct UpdateFunctor {
  UpdateFunctor(DataArray d, double dt) : data(d), dt(dt) {};
  
  KOKKOS_INLINE_FUNCTION
  void operator()(const uint i) const {
    // Mise à jour des vitesses
    data(i, IVX) += dt*data(i, IAX);
    data(i, IVY) += dt*data(i, IAY);
    data(i, IVZ) += dt*data(i, IAZ);

    // Mise à jour des positions
    data(i, IX) += dt * data(i, IVX);
    data(i, IY) += dt * data(i, IVY);
    data(i, IZ) += dt * data(i, IVZ);  
  }

  static void apply(DataArray d, uint N, double dt) {
    UpdateFunctor functor(d, dt);
    Kokkos::RangePolicy policy(0, N);
    Kokkos::parallel_for(policy, functor);
  }

  // Public members
  DataArray data;  
  double dt;
};




Particles::Particles(std::vector<PartHost> &input, 
                     const Params &params) 
    : N(input.size()), params(params) {
  std::cout << "Creating data with " << N << " particles" << std::endl;
  data   = DataArray("Particles", N, FIELD_COUNT);

  double unit_L=1.0, unit_V=1.0, unit_M=1.0;
  if (params.convert_input) {
    unit_L = params.unit_L;
    unit_V = params.unit_V;
    unit_M = params.unit_M;
  }

  for (uint i=0; i < N; ++i) {
    data(i, IX)  = input[i][IX]  / unit_L;
    data(i, IY)  = input[i][IY]  / unit_L;
    data(i, IZ)  = input[i][IZ]  / unit_L;
    data(i, IVX) = input[i][IVX] / unit_V;
    data(i, IVY) = input[i][IVY] / unit_V;
    data(i, IVZ) = input[i][IVZ] / unit_V;
    data(i, IM)  = input[i][IM]  / unit_M;
  }
}

void Particles::ResetAccelerations() {
  ResetAccelerationsFunctor::apply(data, N);
}

void Particles::ComputeAccelerations(const Cell &tree) {
  ComputeAccelerationsFunctor::apply(data, tree, N, 
                                     params.G, params.theta_crit);
}

void Particles::Update(double dt) {
  UpdateFunctor::apply(data, N, dt);
}

}