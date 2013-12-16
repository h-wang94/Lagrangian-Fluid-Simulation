#ifndef _ParticleSystem_h_
#define _ParticleSystem_h_
#include <vector>
#include <iostream>
#include "Particle.h"
#include "SpatialGrid.h"

class ParticleSystem {
  public:
    ParticleSystem();
	  ParticleSystem(Vector grav);
    ParticleSystem(Vector grav, float min_x, float min_y, float min_z, float max_x, float max_y, float max_z);
    void initialize(float timestep);
	  void update(float timestep);
    
    void setBoundaries(float min_x, float min_y, float min_z, float max_x, float max_y, float max_z);

    void addParticle(Particle& p);
    Particle* getParticle(const unsigned int i);
    std::vector<Particle> getParticles();
    float colorFunction(Particle& p);
    std::vector<Particle> getNeighbors(Particle &p);
    Vector surfaceNormal(Particle& p, unsigned const int& i);

  private:
    Vector grav;
    std::vector<Particle> particles;

    float MIN_X;
    float MIN_Y;
    float MIN_Z;
    float MAX_X;
    float MAX_Y;
    float MAX_Z;

    SpatialGrid grid;
    int numRowBoxes;

	  void setDensities();
	  void computeForces();
    void computePressure();
    Vector gravityForce(Particle& p);
    // replace pressure and viscosity force calculation functions
    Vector press_visc(Particle&p, unsigned const int& i);

    Vector pressureForce(Particle& p, unsigned const int& i);
    Vector viscosityForce(Particle& p, unsigned const int& i);
    Vector tensionForce(Particle& p, unsigned const int& i);
    float curvature(Particle& p, unsigned const int& i);

    float defaultKernel(Vector r, const float h);
    Vector gradientKernel(Vector r, const float h);
    float laplacianKernel(Vector r, const float h);
    Vector pressGradientKernel(Vector r, const float h);
    float viscLaplacianKernel(Vector r, const float h);

    void initializeLeapFrog(const float& dt);
    void leapFrog(const float& dt);

    void checkBoundary(Particle& p, Point3D* position, Vector* velocity, Vector* velocityHalf);
};

#endif
