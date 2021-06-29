#include <bits/stdc++.h>

#include <Kokkos_Core.hpp>

#include "particles.h"

using namespace gravity;

constexpr uint saveFreq = 50;
constexpr uint nzeros = 6;

// Gravitational constant
constexpr double dt   = 1.0e-2;
constexpr double tmax = 2.5;

std::vector<PartHost> input;

void Init(std::string filename) {
  std::ifstream f_in;
  f_in.open(filename);
  // 1ere ligne : G
  double G;
  f_in >> G;

  // 2e ligne : Np Nt
  int Np, Nt;
  f_in >> Np >> Nt;
  uint N = Np + Nt;

  for (uint i=0; i < N; ++i) {
    PartHost p;
    f_in >> p[IX] >> p[IY] >> p[IZ] >> p[IVX] >> p[IVY] >> p[IVZ] >> p[IM];
    input.push_back(p);
  }
}

void SaveSolution(uint iteration, Particles &p) {
  std::ostringstream oss;
  oss << "gravity_" << std::setfill('0') << std::setw(nzeros) << iteration << ".3D"; // gravity_00000.3D

  std::ofstream f_out;
  f_out.open(oss.str());
  f_out << "X Y Z V" << std::endl;
  
  p.Send2CPU();
  DataArrayHost data = p.data_h;

  for (uint i=0; i < p.N; ++i) {
    f_out << data(i, IX) << " " << data(i, IY) << " " << data(i, IZ) << " 0.0" << std::endl;
  }
  f_out.close();
}

int main(int argc, char **argv) {
  Kokkos::initialize(argc, argv);
  {
    double t = 0.0;
    uint iteration = 0;

    std::cout << " == Initializing system" << std::endl;
    Init(argv[1]);
    Particles p{input};
    input.clear();

    std::cout << " System : " << std::endl;
    std::cout << "   . Massive particles : " << p.N << std::endl;
    
    // Sauvegarde de la condition initiale
    std::cout << " == Saving initial condition" << std::endl;
    SaveSolution(0, p);
    while (t < tmax) {
      iteration++;
      // 0- Remise à zéro des accélérations
      p.ResetAccelerations();

      // 1- Calcul des accélérations
      p.ComputeAccelerations();

      // 2- Leap-frog
      p.Update(dt);

      // 3- Sauvegarde du système
      if (iteration % saveFreq == 0) {
        std::cout << " == Saving at iteration " << iteration << "; t=" << t << std::endl;
        SaveSolution(iteration, p);
      }

      t += dt;
    }
  }
  Kokkos::finalize();
  
  return 0;
}