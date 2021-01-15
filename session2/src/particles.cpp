#include "particles.h"

namespace gravity {

double Particle::G = 6.67e-11;

Particle::Particle(std::initializer_list<double> init) {
  int i=0;
  for (auto v: init) {
    switch (i) {
      case 0: x = v; break;
      case 1: y = v; break;
      case 2: z = v; break;
      case 3: vx = v; break;
      case 4: vy = v; break;
      case 5: vz = v; break;
      case 6: m  = v; break;
    }

    i++;
  }

  isTracer = false;
}

void Particle::ResetAccelerations() {
  ax = 0.0;
  ay = 0.0;
  az = 0.0;
}

void Particle::ComputeAccelerationMassive(Particle &o) {
  double dx, dy, dz;
  // Distance entre les deux particules
  dx = x - o.x;
  dy = y - o.y;
  dz = z - o.z;
  double r_es = sqrt(dx*dx + dy*dy + dz*dz);

  // Gravitation Newtonienne
  double F = G * m * o.m / (r_es*r_es);

  // Principe Fondamental de la Dynamique (2e loi de Newton)
  ax += -F/m * dx/r_es;
  ay += -F/m * dy/r_es;
  az += -F/m * dz/r_es;

  o.ax +=  F/o.m * dx/r_es;
  o.ay +=  F/o.m * dy/r_es;
  o.az +=  F/o.m * dz/r_es; 
}

void Particle::ComputeAccelerationTracer(Particle &o) {
  double dx, dy, dz;
  // Distance entre les deux particules
  dx = x - o.x;
  dy = y - o.y;
  dz = z - o.z;
  double r_es = sqrt(dx*dx + dy*dy + dz*dz);

  // Gravitation Newtonienne
  double F = G * o.m / (r_es*r_es);

  // Principe Fondamental de la Dynamique (2e loi de Newton)
  ax += -F * dx/r_es;
  ay += -F * dy/r_es;
  az += -F * dz/r_es;
}

void Particle::ComputeAcceleration(Particle &o) {
  if (isTracer)
    ComputeAccelerationTracer(o);
  else
    ComputeAccelerationMassive(o);
}

void Particle::Update(double dt) {
  // Leap-frog
  // 1- Mise à jour des vitesses
  vx += dt * ax;
  vy += dt * ay;
  vz += dt * az;

  // 2- Mise à jour des positions
  x += dt * vx;
  y += dt * vy;
  z += dt * vz;
}

}