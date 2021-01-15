#include <bits/stdc++.h>

#include "particles.h"

using namespace gravity;

constexpr uint saveFreq = 100;
constexpr uint nzeros = 6;

std::vector<Particle> particles; // Particules avec masse
std::vector<Particle> tracers;   // Particules sans masse
uint Np;
uint Nt;

constexpr double dt   = 5.0e-2;
constexpr double tmax = 5.0e3;

void Init(std::string filename) {
  std::ifstream f_in;
  f_in.open(filename);
  // 1ere ligne : G
  f_in >> Particle::G;

  // 2e ligne : Np Nt
  f_in >> Np >> Nt;

  // Np lignes : particules avec masse
  for (uint i=0; i < Np; ++i) {
    Particle p;
    f_in >> p.x >> p.y >> p.z >> p.vx >> p.vy >> p.vz >> p.m;
    particles.push_back(p);
  }

  // Nt lignes : tracers
  for (uint i=0; i < Nt; ++i) {
    Particle t;
    f_in >> t.x >> t.y >> t.z >> t.vx >> t.vy >> t.vz;
    t.isTracer = true;
    tracers.push_back(t);
  }

  f_in.close();
}

void SaveSolution(uint iteration) {
  std::ostringstream oss;
  oss << "gravity_" << std::setfill('0') << std::setw(nzeros) << iteration << ".3D"; // gravity_00000.3D

  std::ofstream f_out;
  f_out.open(oss.str());
  f_out << "X Y Z V" << std::endl;
  for (auto &pi: particles)
    f_out << pi.x << " " << pi.y << " " << pi.z << " 0.0" << std::endl;
  for (auto &ti: tracers) 
    f_out << ti.x << " " << ti.y << " " << ti.z << " 0.0" << std::endl;
  f_out.close();
}

int main(int argc, char **argv) {
  double t = 0.0;
  uint iteration = 0;

  std::cout << " == Initializing system" << std::endl;
  Init(argv[1]);
  std::cout << " System : " << std::endl;
  std::cout << "   . Massive particles : " << Np << " " << particles.size() << std::endl;
  std::cout << "   . Tracers : " << Nt << " " << tracers.size() << std::endl;
  
  // Sauvegarde de la condition initiale
  std::cout << " == Saving initial condition" << std::endl;
  SaveSolution(0);

  while (t < tmax) {
    iteration++;

    // 0- Remise à zéro des accélérations
    for (auto &pi: particles)
      pi.ResetAccelerations();
    for (auto &ti: tracers)
      ti.ResetAccelerations();


    // 1- Calcul des accélérations
    for (int i=0; i < Np-1; ++i) {
      Particle &pi = particles[i];
      for (int j=i+1; j < Np; ++j) {
        Particle &pj = particles[j];
        pi.ComputeAcceleration(pj);
      }
    }

    for (int i=0; i < Nt; ++i) {
      Particle &ti = tracers[i];
      for (int j=0; j < Np; ++j) {
        Particle &pj = particles[j];
        ti.ComputeAcceleration(pj);
      }
    }

    // 2- Leap-frog
    for (auto &p: particles)
      p.Update(dt);
    for (auto &ti: tracers)
      ti.Update(dt);

    // 3- Sauvegarde du système
    if (iteration % saveFreq == 0) {
      std::cout << " == Saving at iteration " << iteration << "; t=" << t << std::endl;
      SaveSolution(iteration);
    }

    t += dt;
  }
  
  return 0;
}