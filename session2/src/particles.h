#pragma once

#include <bits/stdc++.h>

namespace gravity {

struct Particle {
  static double G;
  double x, y, z;
  double vx, vy, vz;
  double ax, ay, az;
  double m;
  bool isTracer;

  Particle() = default;
  Particle(std::initializer_list<double> init);

  void ResetAccelerations();

  void ComputeAcceleration(Particle &o);
  void ComputeAccelerationTracer(Particle &o);
  void ComputeAccelerationMassive(Particle &o);

  void Update(double dt);
};


}