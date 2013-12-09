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
    ParticleSystem(Vector grav, float min_x, float min_y, float min_z, float max_x, float max_y, float max_z, float restCoeff);

    void initialize(float timestep);
    void update(float timestep);
    void setBoundaries(float min_x, float min_y, float min_z, float max_x, float max_y, float max_z);
    void setRestCoeff(float restCoeff);

    void addParticle(Particle& p);
    Particle* getParticle(const unsigned int i);
    std::vector<Particle> getParticles();

  private:
    Vector grav;
    float h;
    float hSq;
    std::vector<Particle> particles;
    
    // boundaries
    float MIN_X;
    float MIN_Y;
    float MIN_Z;
    float MAX_X;
    float MAX_Y;
    float MAX_Z;

    float REST_COEFF;

    bool debug;
    SpatialGrid grid;
    int numRowBoxes;

    void setDensities();
    void computeForces();
    void computePressure();
    Vector gravityForce(Particle& p);
    Vector pressureForce(Particle& p, unsigned const int& i);
    Vector viscosityForce(Particle& p, unsigned const int& i);
    Vector tensionForce(Particle& p, unsigned const int& i);
    float colorFunction(Particle& p, unsigned const int& i);
    Vector surfaceNormal(Particle& p, unsigned const int& i);
    float curvature(Particle& p, unsigned const int& i);

    float defaultKernel(Vector r);
    Vector gradientKernel(Vector r);
    float laplacianKernel(Vector r);
    Vector pressGradientKernel(Vector r);
    float viscLaplacianKernel(Vector r);

    void initializeLeapFrog(const float& dt);
    void leapFrog(const float& dt);

    void checkBoundary(Point3D* position, Vector* velocity, Vector* velocityHalf);
    void bouncebackVelocity(Vector* velocity, Vector normal);

};

#endif
