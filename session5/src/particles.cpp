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
  ComputeAccelerationsFunctor(DataArray d, const Cell &tree, double G) : data(d), tree(tree), G(G) {};

  KOKKOS_INLINE_FUNCTION
  void operator()(const uint i) const {
    Vec3 xi;
    xi[IX] = data(i, IX);
    xi[IY] = data(i, IY);
    xi[IZ] = data(i, IZ);

    Vec3 ai = tree.CalculateAcceleration(xi, G);
    data(i, IAX) = ai[IX];
    data(i, IAY) = ai[IY];
    data(i, IAZ) = ai[IZ];
  }

  static void apply(DataArray d, const Cell &tree, uint N, double G) {
    ComputeAccelerationsFunctor functor(d, tree, G);
    Kokkos::RangePolicy policy(0, N);
    Kokkos::parallel_for(policy, functor);
  }

  // Public members
  DataArray   data;  
  const Cell &tree;
  double G;
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




Particles::Particles(std::vector<PartHost> &input) : N(input.size()) {
  std::cout << "Creating data with " << N << " particles" << std::endl;
  data   = DataArray("Particles", N, FIELD_COUNT);

  for (uint i=0; i < N; ++i) {
    data(i, IX)  = input[i][IX];
    data(i, IY)  = input[i][IY];
    data(i, IZ)  = input[i][IZ];
    data(i, IVX) = input[i][IVX];
    data(i, IVY) = input[i][IVY];
    data(i, IVZ) = input[i][IVZ];
    data(i, IM)  = input[i][IM];
  }
}

void Particles::ResetAccelerations() {
  ResetAccelerationsFunctor::apply(data, N);
}

void Particles::ComputeAccelerations(const Cell &tree) {
  ComputeAccelerationsFunctor::apply(data, tree, N, G);
}

void Particles::Update(double dt) {
  UpdateFunctor::apply(data, N, dt);
}

}