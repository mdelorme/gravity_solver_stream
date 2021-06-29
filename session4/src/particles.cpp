#include "particles.h"

namespace gravity {

/**
 * Foncteur 1: Mise à zero des accélérations
 **/
struct ResetAccelerationsFunctor {
  ResetAccelerationsFunctor(DataArray d) : data(d) {};
  
  KOKKOS_INLINE_FUNCTION
  void operator()(const uint i) const {
    printf("%ld %ld %ld\n", data(i, IAX), data(i, IAY), data(i, IAZ));
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
class ComputeAccelerationsFunctorTeams {
private:
  uint nbTeams;

public:
  using index_t       = uint;
  using team_policy_t = Kokkos::TeamPolicy<Kokkos::IndexType<index_t>>;
  using thread_t      = team_policy_t::member_type;

  void setNbTeams(uint nbTeams_) {nbTeams = nbTeams_;}; 

  ComputeAccelerationsFunctorTeams(DataArray d, uint N) : data(d), N(N) {};

  static void apply(DataArray d, uint N) {
    // Building the functor
    ComputeAccelerationsFunctorTeams functor(d, N);

    // Setting up execution policy
    uint32_t nbTeams_ = N;
    functor.setNbTeams(nbTeams_);

    team_policy_t policy(nbTeams_,
                         Kokkos::AUTO());
    
    Kokkos::parallel_for(policy, functor);
  }

  KOKKOS_INLINE_FUNCTION
  void function(thread_t member) const {
    uint i = member.league_rank();

    // Internal parallel loop
    uint32_t nbElements = N;
    Kokkos::View<double, 3> acc();
    acc(IX) = 0.0;
    acc(IY) = 0.0;
    acc(IZ) = 0.0;

    Kokkos::parallel_reduce(
      Kokkos::TeamThreadRange(member, nbElements),
      KOKKOS_LAMBDA(const index_t j, auto &local_acc) {
        if (i==j)
          return;

        double dx, dy, dz;
        // Distance entre les deux particules
        dx = data(i, IX) - data(j, IX);
        dy = data(i, IY) - data(j, IY);
        dz = data(i, IZ) - data(j, IZ);
        double r_es = sqrt(dx*dx + dy*dy + dz*dz);

        // Gravitation Newtonienne
        double F = 0.044955 * data(j, IM) / (r_es*r_es);

        // Principe Fondamental de la Dynamique (2e loi de Newton)
        local_acc(IX) -= -F*dx/r_es;
        local_acc(IY) -= -F*dy/r_es;       
        local_acc(IZ) -= -F*dz/r_es;
      }, acc);

    member.team_barrier();

    Kokkos::single(Kokkos::PerTeam(member), [=]() {
      data(IAX, i) = acc(IX);
      data(IAY, i) = acc(IY);
      data(IAZ, i) = acc(IZ);
    });    
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(thread_t member) const {
    function(member);
  }

  // Attributes here
  DataArray data; //!< Données en mémoire GPU
  uint      N;    //!< Nombre de particules
};

struct ComputeAccelerationsFunctorRange {
  ComputeAccelerationsFunctorRange(DataArray d, uint N) : data(d), N(N) {};
  
  KOKKOS_INLINE_FUNCTION
  void operator()(const uint i) const {
    for (uint j=0; j < N; ++j) {
      if (i==j)
        continue;

      double dx, dy, dz;
      // Distance entre les deux particules
      dx = data(i, IX) - data(j, IX);
      dy = data(i, IY) - data(j, IY);
      dz = data(i, IZ) - data(j, IZ);
      double r_es = sqrt(dx*dx + dy*dy + dz*dz);

      // Gravitation Newtonienne
      double F = 0.044955 * data(j, IM) / (r_es*r_es);

      // Principe Fondamental de la Dynamique (2e loi de Newton)
      data(i, IAX) -= F * dx/r_es;
      data(i, IAY) -= F * dy/r_es;
      data(i, IAZ) -= F * dz/r_es;
    }  
  }

  static void apply(DataArray d, uint N) {
    ComputeAccelerationsFunctorRange functor(d, N);
    Kokkos::RangePolicy policy(0, N);
    Kokkos::parallel_for(policy, functor);
  }

  // Public members
  DataArray data;  
  uint N;
};

using ComputeAccelerationsFunctor = ComputeAccelerationsFunctorTeams;

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
    data(i, IX) += dt*data(i, IVX);
    data(i, IY) += dt*data(i, IVY);
    data(i, IZ) += dt*data(i, IVZ);  
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
  data_h = Kokkos::create_mirror(data);

  for (uint i=0; i < N; ++i) {
    data_h(i, IX)  = input[i][IX];
    data_h(i, IY)  = input[i][IY];
    data_h(i, IZ)  = input[i][IZ];
    data_h(i, IVX) = input[i][IVX];
    data_h(i, IVY) = input[i][IVY];
    data_h(i, IVZ) = input[i][IVZ];
    data_h(i, IM)  = input[i][IM];
  }

  Send2GPU();
}

void Particles::ResetAccelerations() {
  ResetAccelerationsFunctor::apply(data, N);
}

void Particles::ComputeAccelerations() {
  ComputeAccelerationsFunctor::apply(data, N);
}

void Particles::Update(double dt) {
  UpdateFunctor::apply(data, N, dt);
}

void Particles::Send2GPU() {
  Kokkos::deep_copy(data, data_h);
}

void Particles::Send2CPU() {
  Kokkos::deep_copy(data_h, data);
}

}