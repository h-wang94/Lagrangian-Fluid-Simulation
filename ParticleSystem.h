#ifndef _ParticleSystem_h_
#define _ParticleSystem_h_
#include <vector>
#include "Particle.h"

class ParticleSystem {
  public:
  
	  ParticleSystem(float grav);
	  void computeForces();
    void computePressure();
	  void update(float timestep);
	  void setDensities();
    

  private:
    std::vector<Particle> particles;

    Vector gravityForce(Particle& p);
    Vector pressureForce(Particle& p);
    Vector viscosityForce(Particle& p);
    float defaultKernel(Vector r, float h);
    /*Vector gradientKernel(Vector r, float h);
    float laplacianKernel(Vector r, float h);*/
    Vector pressGradientKernel(Vector r, float h);
    float viscLaplacianKernel(Vector r, float h);
    void leapFrog(float dt);
    Vector grav;
};

#endif
