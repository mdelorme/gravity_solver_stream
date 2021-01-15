#include <bits/stdc++.h>

constexpr double G  = 39.147;       // AU^3 / Msun / yr^2
constexpr double Ms = 1.0;          // Msun
constexpr double Me = 3.0025e-6;    // Msun
constexpr double AU = 1.0;          // Au
constexpr double dt = 1.0 / 365.25; // 1j en yr

// Terre
double xe, ye, ze; 
double vxe, vye, vze;
double axe, aye, aze;

// Soleil
double xs, ys, zs;
double vxs, vys, vzs;
double axs, ays, azs;

int main(int argc, char **argv) {
  // Initialisation
  // Soleil
  xs  = 0.0;
  ys  = 0.0;
  zs  = 0.0;
  vxs = 0.0;
  vys = 0.0;
  vzs = 0.0;
  
  // Terre
  xe = AU;
  ye = 0.0;
  ze = 0.0;
  vxe = 0.0;
  vye = sqrt(G * Ms / xe);
  vze = 0.0;

  double t = 0.0;
  double tmax = 10.0;

  while (t < tmax) {
    // 1- Calcul des accelerations
    double dx, dy, dz;
    // vecteur S -> E
    dx = xe - xs;
    dy = ye - ys;
    dz = ze - zs;
    double r_es = sqrt(dx*dx + dy*dy + dz*dz);

    // Gravitation Newtonienne
    double F = G * Me * Ms / (r_es*r_es);

    // Principe Fondamental de la Dynamique (2e loi de Newton)
    axe = -F/Me * dx/r_es;
    aye = -F/Me * dy/r_es;
    aze = -F/Me * dz/r_es;

    axs =  F/Ms * dx/r_es;
    ays =  F/Ms * dy/r_es;
    azs =  F/Ms * dz/r_es; 

    // 2- Leap-frog
    // a- Mise à jour des vélocités
    vxe = vxe + dt * axe;
    vye = vye + dt * aye;
    vze = vze + dt * aze;

    vxs = vxs + dt * axs;
    vys = vys + dt * ays;
    vzs = vzs + dt * azs; 

    // b- Mise à jour des positions
    xe = xe + dt * vxe;
    ye = ye + dt * vye;
    ze = ze + dt * vze;

    xs = xs + dt * vxs;
    ys = ys + dt * vys;
    zs = zs + dt * vzs;
    

    // 3- Affichage de la position de la Terre
    std::cout << t << " " << xe << " " << ye << std::endl;

    t += dt;
  }
  
  return 0;
}