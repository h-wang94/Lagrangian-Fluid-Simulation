#ifndef _ParticleSystem_h_
#define _ParticleSystem_h_
#include <vector>
#include "Particle.h"

class ParticleSystem {
  public:
  
	  ParticleSystem(float grav);
	  void computeForces();

  private:
    std::vector<Particle> particles;

};

#endif
