#pragma once

#include <bits/stdc++.h>
#include "params.h"

namespace gravity {


KOKKOS_INLINE_FUNCTION
double dist(Vec3 v1, Vec3 v2) {
  double d = 0.0;
  for (int dir=0; dir < 3; ++dir)
    d += (v1[dir]-v2[dir])*(v1[dir]-v2[dir]);
  return sqrt(d);
}

KOKKOS_INLINE_FUNCTION
Vec3 operator+(Vec3 v1, Vec3 v2) {
  return Vec3{v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2]};
}

KOKKOS_INLINE_FUNCTION
Vec3 operator+(Vec3 v1, double q) {
  return Vec3{v1[0]+q, v1[1]+q, v1[2]+q};
}

KOKKOS_INLINE_FUNCTION
Vec3 operator+=(Vec3 &v1, const Vec3 &v2) {
  v1[0] += v2[0];
  v1[1] += v2[1];
  v1[2] += v2[2];

  return v1;
}

/**
 * 1 cellule de l'octree
 **/
struct Cell {
  Cell() = default;
  Cell(Vec3 x0, double size);
  ~Cell() = default;

  void DistributeParticles(const std::vector<TreePart> &particles, int level=0);
  int PrintStats(int level=0) const;

  KOKKOS_INLINE_FUNCTION
  double CalculateOpeningAngle(const Vec3 xi) const {
    return size / dist(xi, x0+hsize);
  }

  KOKKOS_INLINE_FUNCTION
  Vec3 CalculateAcceleration(const Vec3 xi, const double G, const double theta_crit) const {
    Vec3 ai {0.0, 0.0, 0.0};

    double theta = CalculateOpeningAngle(xi);
    if (theta < theta_crit || Np == 1) {
      Vec3 acell;
      const double d = dist(xi, cm);
      if (d < 1.0e-10)
        return ai;

      const double d3 = d*d*d; 
      double F = G * total_mass / d3;
      ai[IX] += F * (cm[IX] - xi[IX]);
      ai[IY] += F * (cm[IY] - xi[IY]);
      ai[IZ] += F * (cm[IZ] - xi[IZ]);
    }
    else if (children.size() > 0) {
      for (int i=0; i < 8; ++i)
        ai += children[i].CalculateAcceleration(xi, G, theta_crit);
    }

    return ai;
  }

  double total_mass; // !< Masse totale contenue dans la cellule
  Vec3   cm;         // !< Centre de masse
  double size;       // !< Taille de la cellule
  double hsize;      // !< Taille d'une demi-cellule
  int    Np;         // !< Nombre de particules
  Vec3   x0;         // !< Bounding-Box : x0 = coin inférieur-gauche

  std::vector<Cell> children; // !< 8 enfants "if any"
};

}