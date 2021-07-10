#include "tree.h"

namespace gravity {

Cell::Cell(Vec3 x0, double size)
  : x0(x0), size(size), hsize(size*0.5) {};

void Cell::DistributeParticles(const std::vector<TreePart> &particles, int level) {
  // 3 cas :

  // 1- particles.size() = 0 -> sort rien à faire
  // 2- particles.size() = 1 -> on définit le centre de masse, et la masse tot et on sort
  // 3- particles.size() > 1 -> on coupe en 8, on calcule les vecteur<TreePart> pour chaque quadrant
  //                            + cdm + masse tot, appels récursifs sur DistributeParticles

  Np = particles.size();

  if (Np == 1) {
    cm[IX] = particles[0][IX];
    cm[IY] = particles[0][IY];
    cm[IZ] = particles[0][IZ];
    total_mass = particles[0][IM];
  }
  else if (Np > 1) {
    cm[IX] = 0.0;
    cm[IY] = 0.0;
    cm[IZ] = 0.0;
    total_mass = 0.0;

    // Création des 8 enfants de ma cellule courante
    for (int i=0; i < 8; ++i) {
      Vec3 xi0;
      int ii = i/4;               
      int ij = (i-ii*4)/2;
      int ik = i-ii*4-ij*2;
      xi0[IX] = x0[IX] + ii*hsize;
      xi0[IY] = x0[IY] + ij*hsize;
      xi0[IZ] = x0[IZ] + ik*hsize;
      children.push_back(Cell{xi0, hsize});
    }

    // On remplit les sous-vecteurs de particules
    // Et on calcule le centre de masse et la masse totale
    std::vector<TreePart> sub_particles[8];
    for (auto &p: particles) {
      // Définition du centre de masse et de la masse totale
      cm[IX] += p[IX] * p[IM];
      cm[IY] += p[IY] * p[IM];
      cm[IZ] += p[IZ] * p[IM];
      total_mass += p[IM];

      int ii = (p[IX] - x0[IX]) / hsize;
      int ij = (p[IY] - x0[IY]) / hsize;
      int ik = (p[IZ] - x0[IZ]) / hsize;

      int cid = ii*4 + ij*2 + ik;
      sub_particles[cid].push_back(p);
    }

    // On divise le cm par la masse totale
    if (total_mass > 0.0)
      for (int dir=0; dir < 3; ++dir)
        cm[dir] /= total_mass;

    // Appel récursif à Distribute Particles
    for (int i=0; i < 8; ++i)
      children[i].DistributeParticles(sub_particles[i], level+1);
  }
  else {
    cm[IX] = 0.0;
    cm[IY] = 0.0;
    cm[IZ] = 0.0;
    total_mass = 0.0;
  }
}

std::map<int, int> cells_per_level;

int Cell::PrintStats(int level) const {
  if (level == 0)
    cells_per_level.clear();

  if (cells_per_level.count(level) == 0)
    cells_per_level[level] = 1;
  else
    cells_per_level[level]++;

  int max_level = level;
  if (children.size() > 0) {
    for (auto &c : children) {
      int max_lev = c.PrintStats(level+1);
      max_level = std::max(max_level, max_lev);
    }
  }
  if (level == 0) {
    std::cout << "==== Cells per level ====" << std::endl;
    for (int i=0; i <= max_level; ++i)
      std::cout << " Level " << i << "\thas " << cells_per_level[i] << std::endl;
  }

  return max_level;
}

}