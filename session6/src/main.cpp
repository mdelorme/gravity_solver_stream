#include <bits/stdc++.h>

#include <Kokkos_Core.hpp>

#include "particles.h"
#include "tree.h"
#include "params.h"

using namespace gravity;


// Vecteur de particules lues depuis la ligne de commande
std::vector<PartHost> input;
// Arbre de Barnes & Hut
Cell tree;

Params Init(std::string filename) {
  Params params{filename};

  std::ifstream f_in;
  std::string input_filename = "../datasets/" + params.input_filename;
  f_in.open(input_filename);
  std::cout << input_filename << std::endl;

  while(f_in.good()) {
    PartHost p;
    f_in  >> p[IX] >> p[IY] >> p[IZ] 
          >> p[IVX] >> p[IVY] >> p[IVZ]
          >> p[IM];     
    if (!f_in.good())
      break;

    input.push_back(p);
  }

  return params;
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

void SaveSolution(uint iteration, Particles &p, const Params &params) {
  std::ostringstream oss;
  oss << params.output_prefix + "_" << std::setfill('0') << std::setw(params.nzeros) << iteration << ".3D"; // output_prefix_00000.3D

  std::ofstream f_out;
  f_out.open(oss.str());
  f_out << "X Y Z V" << std::endl;

  for (uint i=0; i < p.N; ++i) {
    f_out << p.data(i, IX) * params.unit_L << " " 
          << p.data(i, IY) * params.unit_L << " " 
          << p.data(i, IZ) * params.unit_L << " 0.0" << std::endl;
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
    Params params = Init(argv[1]);
    Particles p{input, params};
    input.clear();

    std::cout << " System : " << std::endl;
    std::cout << "   . Massive particles : " << p.N << std::endl;
    
    // Sauvegarde de la condition initiale
    std::cout << " == Saving initial condition" << std::endl;
    SaveSolution(0, p, params);
    while (t < params.tmax) {
      iteration++;
      // 0- Remise à zéro des accélérations
      p.ResetAccelerations();

      // 1- Construction de l'arbre
      BuildTree(p);

      // 2- Calcul des accélérations
      p.ComputeAccelerations(tree);

      // 3- Leap-frog
      p.Update(params.dt);

      // 4- Sauvegarde du système
      if (iteration % params.saveFreq == 0) {
        std::cout << " == Saving at iteration " << iteration << "; t=" << t << std::endl;
        SaveSolution(saveId, p, params);
        saveId++;
        //tree.PrintStats();
      }

      t += params.dt;
    }
  }
  Kokkos::finalize();
  
  return 0;
}