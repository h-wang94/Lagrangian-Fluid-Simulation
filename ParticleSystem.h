#ifndef _ParticleSystem_h_
#define _ParticleSystem_h_
#include <vector>
#include <iostream>
#include "Particle.h"

class ParticleSystem {
  public:
    ParticleSystem();
	  ParticleSystem(Vector grav);

    void initialize(float timestep);
	  void update(float timestep);

    void addParticle(Particle& p);
    Particle* getParticle(const unsigned int i);
    std::vector<Particle> getParticles();

  private:
    Vector grav;
    float h;
    float hSq;
    bool debug;

	  void setDensities();
	  void computeForces();
    void computePressure();
    Vector gravityForce(Particle& p);
    Vector pressureForce(Particle& p, unsigned const int& i);
    Vector viscosityForce(Particle& p, unsigned const int& i);

    float defaultKernel(Vector r);
    Vector pressGradientKernel(Vector r);
    float viscLaplacianKernel(Vector r);

    void initializeLeapFrog(const float& dt);
    void leapFrog(const float& dt);

    void checkBoundary(Point3D* position, Vector* velocity, Vector* velocityHalf);

    std::vector<Particle> particles;
};

#endif
