#include <bits/stdc++.h>

#include <Kokkos_Core.hpp>

#include "particles.h"
#include "tree.h"
#include "params.h"

using namespace gravity;

constexpr uint saveFreq = 100;
constexpr uint nzeros = 6;

// Gravitational constant
constexpr double dt   = 1.0e-1;
constexpr double tmax = 10000.0;

std::vector<PartHost> input;
Cell tree;
double G;

void Init(std::string filename) {
  std::ifstream f_in;
  f_in.open(filename);
  // 1ere ligne : G
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

void BuildTree(const Particles &p) {
  // 1- Construction de la liste des particules et calcul de la bounding box
  Vec3 x0, x1;
  double size = 0.0;
  std::vector<TreePart> particles;

  for (int i=0; i < p.N; ++i) {
    TreePart pi;
    pi[IX] = p.data(i, IX);
    pi[IY] = p.data(i, IY);
    pi[IZ] = p.data(i, IZ);
    pi[IM] = p.data(i, IM);

    particles.push_back(pi);

    // Calcul de la bounding-box
    if (i==0) {
      x0[IX] = pi[IX];
      x1[IX] = pi[IX];
      x0[IY] = pi[IY];
      x1[IY] = pi[IY];
      x0[IZ] = pi[IZ];
      x1[IZ] = pi[IZ];
    }
    else {
      for (int dir=0; dir < 3; ++dir) {
        if (pi[dir] < x0[dir]) {
          x0[dir] = pi[dir];
          size = std::max(size, x1[dir]-x0[dir]);
        }
        else if (pi[dir] > x1[dir]) {
          x1[dir] = pi[dir];
          size = std::max(size, x1[dir]-x0[dir]);
        }
      }
    }
  }

  double epsilon = 1.0e-3;
  for (int i=0; i < 3; ++i)
    x0[i] -= epsilon;
  size += 2*epsilon;

  // 2- Création de la "Root"
  tree = Cell(x0, size);
  tree.DistributeParticles(particles);
}

void SaveSolution(uint iteration, Particles &p) {
  std::ostringstream oss;
  oss << "gravity_" << std::setfill('0') << std::setw(nzeros) << iteration << ".3D"; // gravity_00000.3D

  std::ofstream f_out;
  f_out.open(oss.str());
  f_out << "X Y Z V" << std::endl;

  for (uint i=0; i < p.N; ++i) {
    f_out << p.data(i, IX) << " " << p.data(i, IY) << " " << p.data(i, IZ) << " 0.0" << std::endl;
  }
  f_out.close();
}

int main(int argc, char **argv) {
  Kokkos::initialize(argc, argv);
  {
    double t = 0.0;
    uint iteration = 0;
    uint saveId = 0;

    std::cout << " == Initializing system" << std::endl;
    Init(argv[1]);
    Particles p{input};
    p.G = G;
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

      // 1- Construction de l'arbre
      BuildTree(p);

      // 2- Calcul des accélérations
      p.ComputeAccelerations(tree);

      // 3- Leap-frog
      p.Update(dt);

      // 4- Sauvegarde du système
      if (iteration % saveFreq == 0) {
        std::cout << " == Saving at iteration " << iteration << "; t=" << t << std::endl;
        SaveSolution(saveId, p);
        saveId++;
        //tree.PrintStats();
      }

      t += dt;
    }
  }
  Kokkos::finalize();
  
  return 0;
}