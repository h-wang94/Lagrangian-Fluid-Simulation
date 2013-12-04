#include "Particle.h"
#include "Vector.h"
#include "Point3D.h"
#include "ParticleSystem.h"
#include <iostream>
#include <cstdlib>
using namespace std;
int main() {
  cout << "//===========================================//" << endl
       << "// Fluid Simulation Project " << endl
       << "//===========================================//" << endl;

  // Testing print/debug statement
  /*Particle particle = Particle(1, 2, 3, Vector(4, 5, 6), Vector(7, 8, 9));*/
  /*cout << particle;*/
  ParticleSystem system = ParticleSystem(Vector(0,0,-9.8));
  Particle p1 = Particle(1, 1, 1, 1, Vector(0,0,1), Point3D(5, 5, 5));
  //Particle p2 = Particle(0.02, 998.29, 1, 1, Vector(0,1,1), Point3D(1, 1, 1));
  for(unsigned int i = 0; i < 2; i++) {
    system.addParticle(p1);
    //system.addParticle(p2);
  }
  system.initialize();
  for(float i = 0; i < 0.1; i+=0.1) {
    system.update(i);
  }
}
