#pragma once

#include <Kokkos_Core.hpp>
#include <bits/stdc++.h>

#include "../../external/inih/INIReader.h"

namespace gravity {

constexpr double G_cgs = 6.67e-8;

struct Params {
  //!< Critère d'ouverture pour l'arbre
  double theta_crit;

  //!< Fréquence de sauvegarde en itérations
  uint saveFreq;

  //!< Nombre de zéros de padding dans le nom des sauvegardes
  uint nzeros;

  //!< Temps maximum de simulation
  double tmax = 1.0;

  //!< Pas de temps
  double dt;

  //!< Constante gravitationnelle
  double G;

  //!< Unité de longueur
  double unit_L;

  //!< Unité de masse
  double unit_M;

  //!< Unité de temps
  double unit_T;

  //!< Unité de vitesse
  double unit_V;

  //!< Est-ce que le fichier d'input est déjà normalisé
  bool convert_input;

  //!< Fichier d'input
  std::string input_filename;

  //!< Prefixe output
  std::string output_prefix;

  Params(std::string filename) {
    INIReader reader(filename);

    tmax = reader.GetReal("run", "tmax", 1.0);
    dt   = reader.GetReal("run", "dt",   1.0e-2);
    input_filename = reader.Get("run", "input_filename", "");
    if (input_filename == "") {
      // ERROR
    }

    saveFreq = reader.GetInteger("saves", "saveFreq", 100);
    nzeros   = reader.GetInteger("saves", "nzeros", 4);
    output_prefix = reader.Get("saves", "output_prefix", input_filename);

    theta_crit = reader.GetReal("barnes-hut", "theta_crit", 0.1);

    // Lecture des unités
    G      = reader.GetReal("physics", "G", 1.0);
    unit_L = reader.GetReal("physics", "length_unit", 1.0);
    unit_M = reader.GetReal("physics", "mass_unit",   1.0);
    unit_T = sqrt(G / G_cgs * std::pow(unit_L, 3.0) / unit_M);
    unit_V = unit_L / unit_T;
    convert_input = reader.GetBoolean("physics", "convert_input", true);
  }
};

// Particle table : double[N][10]
enum Fields : uint8_t {
  IX  = 0,
  IY  = 1,
  IZ  = 2,
  IM  = 3,
  IVX = 4,
  IVY = 5,
  IVZ = 6,
  IAX = 7,
  IAY = 8,
  IAZ = 9,

  FIELD_COUNT=10
};

using DataArray = Kokkos::View<double**>;
using PartHost  = std::array<double, FIELD_COUNT>;
using TreePart  = std::array<double, 4>;
using Vec3      = std::array<double, 3>;

}