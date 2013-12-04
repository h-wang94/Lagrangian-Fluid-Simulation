#ifndef _ParticleSystem_h_
#define _ParticleSystem_h_
#include <vector>
#include "Particle.h"

class ParticleSystem {
  public:
  
	  ParticleSystem(Vector grav);
    void initialize();
	  void computeForces();
    void computePressure(const float stiffness, const float restDensity);
	  void update(float timestep);
	  void setDensities();
    void addParticle(Particle& p);
    std::vector<Particle> particles;
    

  private:

    Vector gravityForce(Particle& p);
    Vector pressureForce(Particle& p, unsigned const int i);
    Vector viscosityForce(Particle& p, unsigned const int i);
    float defaultKernel(Vector r, const float h);
    /*Vector gradientKernel(Vector r, float h);
    float laplacianKernel(Vector r, float h);*/
    Vector pressGradientKernel(Vector r, const float h);
    float viscLaplacianKernel(Vector r, const float h);
    void leapFrog(const float dt);
    void initializeLeapFrog(const float dt);
    void checkBoundary(Point3D* position, Vector* velocity);
    Vector grav;
};

#endif
