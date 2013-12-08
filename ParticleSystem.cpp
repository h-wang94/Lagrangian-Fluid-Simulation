#include "ParticleSystem.h"
#include <omp.h>
#include <time.h>

#define PI 3.14159265
#define MAX_X .1
#define MAX_Y 1
#define MAX_Z .1
#define MIN_X -.1
#define MIN_Y -1
#define MIN_Z -.1
#define REST_COEFF 0.8

ParticleSystem::ParticleSystem() {
  this->grav = Vector(0,0,-9.8);
  this->h = 0.0457;
  this->hSq = pow(h, 2.0f);
  this->debug = true;
}

ParticleSystem::ParticleSystem(Vector grav){
  this->grav = grav;
  this->h = 0.0457; // doesnt seem to do much interaction for 100ish particles
  //this->h = 1; // for funky fusion
  this->hSq = pow(h, 2.0f);
  this->debug = true;
  this->numRowBoxes = 100;
  this->grid = SpatialGrid(100, h);
}

void ParticleSystem::initialize(float timestep) {
  this->setDensities();
  this->computePressure();
  this->computeForces();
  this->initializeLeapFrog(-timestep);
}

void ParticleSystem::update(float timestep){
  this->setDensities();//for each particle, compute particle's density
  this->computePressure(); // now compute each particle's pressure. randomly put in numbers.
  this->computeForces();
  this->leapFrog(timestep);
  grid.updateBoxes(particles);
}

Particle* ParticleSystem::getParticle(const unsigned int i) {
  return &particles[i];
}

std::vector<Particle> ParticleSystem::getParticles() {
  return this->particles;
}

// computes the internal and external forces.
// also sets acceleration
void ParticleSystem::computeForces(){
  Vector gravity, pressure, viscosity, tension, force;
#pragma omp parallel for firstprivate(gravity, pressure, viscosity, tension, force)
  for(unsigned int i = 0; i < particles.size(); i++) {
    gravity = gravityForce(particles[i]);
    pressure = pressureForce(particles[i], i);
    viscosity = viscosityForce(particles[i], i);
    tension = tensionForce(particles[i], i);

    force = gravity - pressure + viscosity;
    particles[i].setAcceleration(force / particles[i].getDensity()); // intuitively using mass..but slides say density...
  }
}

void ParticleSystem::computePressure() {
  float pressure = 0;
#pragma omp parallel for firstprivate(pressure)
  for(unsigned int i = 0; i < particles.size(); i++) {
    // p = k ( (p / p0)^7 - 1)
    pressure = particles[i].getStiffness() * (pow((particles[i].getDensity() / particles[i].getRestDensity()), 7.0f) - 1.0f);
    particles[i].setPressure(pressure);
  }
}

// need to look at this
void ParticleSystem::setDensities(){

  float density = 0;
#pragma omp parallel for firstprivate(density)
  for(unsigned int i = 0; i < particles.size(); i++){
    density = 0;

    /*std::vector<Particle> list = grid.getNeighbors(particles[i]);
    for(unsigned int j = 0; j < list.size(); j++){ 
      Vector dist = particles[i].getPosition() - list[j].getPosition();
      //if (dist.getMagnitude() <= hSq) {
      if (dist.getMagnitude() <= h) {
        density += defaultKernel(dist) * list[j].getMass();
      }*/

      for(unsigned int j = 0; j < particles.size(); j++){ // need to use spatial grid	
        Vector dist = particles[i].getPosition() - particles[j].getPosition();
      //if (dist.getMagnitude() <= hSq) {
      if (dist.getMagnitude() <= h) {
      density += defaultKernel(dist) * particles[j].getMass();
      }

    }
    if(density == 0){
      density = particles[i].getMass() / .00000001;
    }

    particles[i].setDensity(density);											//set the particle[i]'s density to particle[i]
  }
}

Vector ParticleSystem::gravityForce(Particle& p) {
  return grav * p.getDensity(); 
}

Vector ParticleSystem::pressureForce(Particle& p, unsigned const int& i) {
  Vector pressure;
  float coeff = 0;
  unsigned int j = 0;
  // this is so stupid...........but ill think of a better way. i dont wanna do checks cause it might make a difference since this is computed every time for every particle
#pragma omp parallel for firstprivate(pressure, coeff)
  for(j = 0; j < i; j++) {
    Vector diff = p.getPosition() - particles[j].getPosition();
    //if (diff.getMagnitude() <= hSq) {
    if (diff.getMagnitude() <= h) {
      coeff = (p.getPressure() + particles[j].getPressure()) / 2.0f * particles[j].getVolume();
      pressure += pressGradientKernel(diff) * coeff;
    }
  }
 #pragma omp parallel for firstprivate(pressure, coeff)
  for(j = j + 1; j < particles.size(); j++) {
    Vector diff = p.getPosition() - particles[j].getPosition();
    //if (diff.getMagnitude() <= hSq) {
    if (diff.getMagnitude() <= h) {
      coeff = (p.getPressure() + particles[j].getPressure()) / 2.0f * particles[j].getVolume();
      pressure += pressGradientKernel(diff) * coeff;
    }
  }
  return pressure;
}

Vector ParticleSystem::viscosityForce(Particle& p, unsigned const int& i) {
  Vector viscosity;
  Vector coeff;
  unsigned int j;
  for(j = 0; j < i; j++) {
    Vector diff = p.getPosition() - particles[j].getPosition();
    //if (diff.getMagnitude() <= hSq) {
    if (diff.getMagnitude() <= h) {
      coeff = (particles[j].getVelocity() - p.getVelocity()) * particles[j].getVolume();
      viscosity += coeff * viscLaplacianKernel(diff);
    }
  }
  for(j = j + 1; j < particles.size(); j++) {
    Vector diff = p.getPosition() - particles[j].getPosition();
    //if (diff.getMagnitude() <= hSq) {
    if (diff.getMagnitude() <= h) {
      coeff = (particles[j].getVelocity() - p.getVelocity()) * particles[j].getVolume();
      viscosity += coeff * viscLaplacianKernel(diff);
    }
  }
  viscosity = viscosity * p.getViscosity(); 
  return viscosity;
  }

  Vector ParticleSystem::tensionForce(Particle& p, unsigned const int& i) {
    Vector normal = surfaceNormal(p, i);
    if (normal.getMagnitude() > 7.065) { // threshold
      normal.normalize();
      return normal * curvature(p, i) * (-0.0728); // surface tension constant for water
    }
    return Vector(0, 0, 0);
  }

float ParticleSystem::colorFunction(Particle& p, unsigned const int& i) {
  float color = 0;
  Vector diff;
  unsigned int j;
  for(j = 0; j < i; j++) {
    diff = p.getPosition() - particles[j].getPosition();
    if (diff.getMagnitude() <= h) {
      color += particles[j].getMass() / particles[j].getDensity() * defaultKernel(diff);
    }
  }
  for(j = j + 1; j < particles.size(); j++) {
    diff = p.getPosition() - particles[j].getPosition();
    if (diff.getMagnitude() <= h) {
      color += particles[j].getMass() / particles[j].getDensity() * defaultKernel(diff);
    }
  }
  return color;
}

Vector ParticleSystem::surfaceNormal(Particle& p, unsigned const int& i) {
  Vector normal = Vector(0, 0, 0);
  Vector diff;
  unsigned int j;
  for(j = 0; j < i; j++) {
    diff = p.getPosition() - particles[j].getPosition();
    if (diff.getMagnitude() <= h) {
      normal += gradientKernel(diff) * particles[j].getMass() / particles[j].getDensity();
    }
  }
  for(j = j + 1; j < particles.size(); j++) {
    diff = p.getPosition() - particles[j].getPosition();
    if (diff.getMagnitude() <= h) {
      normal += gradientKernel(diff) * particles[j].getMass() / particles[j].getDensity();
    }
  }
  return normal;
}

float ParticleSystem::curvature(Particle& p, unsigned const int& i) {
  float curvature = 0;
  Vector diff;
  unsigned int j;
  for(j = 0; j < i; j++) {
    diff = p.getPosition() - particles[j].getPosition();
    if (diff.getMagnitude() <= h) {
      curvature += particles[j].getMass() / particles[j].getDensity() * laplacianKernel(diff);
    }
  }
  for(j = j + 1; j < particles.size(); j++) {
    diff = p.getPosition() - particles[j].getPosition();
    if (diff.getMagnitude() <= h) {
      curvature += particles[j].getMass() / particles[j].getDensity() * laplacianKernel(diff);
    }
  }
  return curvature;
}


// Smoothing distance h is half of the difference between particle i's most distant nearest neighbor and i
// Poly6 Kernels used for everything except pressure and viscosity forces 
// Notes: http://image.diku.dk/projects/media/kelager.06.pdf (Page 16)
// Less expensive compared to gaussian one because of computation of e and no square roots
// Not sure whether I'm supposed to normalize Vector r though

float ParticleSystem::defaultKernel(Vector r) {
  float rMag = r.getMagnitude();
  return (315.0f * pow((pow(h, 2.0f) - pow(rMag, 2.0f)),3.0) / (64.0f * PI * pow(h, 9.0f)));
}

/* gradient and laplacian of poly6 kernels. prob not needed if we use the spiky and viscosity kernels for other calculations*/
Vector ParticleSystem::gradientKernel(Vector r) {
  float rMag = r.getMagnitude();
  float coeff = pow(pow(h, 2.0f) - pow(rMag, 2.0f), 2.0f) * -945 / (32 * PI * pow(h, 9.0f));
  return r * coeff;
}

float ParticleSystem::laplacianKernel(Vector r) {
  float rMag = r.getMagnitude();
  return (-945 / (32 * PI * pow(h, 9.0f))) * (pow(h, 2.0f) - pow(rMag, 2.0f)) * (3 * pow(h, 2.0f) - 7 * pow(rMag, 2.0f));
}

// Spiky Kernel to calculate pressure 
// Give more repulsive pressure at float distance and thus avoids clustering.
Vector ParticleSystem::pressGradientKernel(Vector r) {
  float rMag = r.getMagnitude();
  /*if (rMag == 0) {
    return Vector(0,0,0);
  }*/
  float coeff = (-45 * pow((h - rMag), 2.0f)) / (PI * pow(h, 6.0f));
  //r.normalize();
  return r * coeff;
}

// Viscosity kernel to calculate viscosity
// Gives more stable viscosity and makes it possible to damp simulation better
// Laplacian in poly6 kernel becomes negative really fast. The viscosity kernel's laplacian is positive everywhere
float ParticleSystem::viscLaplacianKernel(Vector r) {
  float rMag = r.getMagnitude();
  return (45 * (h - rMag)) / (PI * pow(h, 6.0f));

}

// Leap frog integration. Takes in a dt and will loop through all the particles in our system.
// Can be changed to work with leapfrogging just a certain particle.
void ParticleSystem::leapFrog(const float& dt) {
  for(unsigned int i = 0; i < particles.size(); i++) {
    Particle& p = particles[i];
    // get velocity at time t - dt/2. v_{t - dt/2}
    Vector old = p.getVelocityHalf(); 
    // velocity at time t + dt/2. v_{t + dt/2} = v_{t - dt / 2} + a * dt
    Vector tempVelocityHalf = old + p.getAcceleration() * dt;
    // set position at time t. pos_{t + dt} = pos_{t} + v_{t + dt / 2} * dt
    Point3D tempPosition = p.getPosition() + (tempVelocityHalf * dt);
    Vector tempVelocity = (old + p.getVelocityHalf()) / 2.0f;
    checkBoundary(&tempPosition, &tempVelocity, &tempVelocityHalf);
    p.setOldPosition(p.getPosition());
    p.setVelocityHalf(tempVelocityHalf); 
    p.setPosition(tempPosition); 
    // use midpoint approximation for velocity at time t. v_{t} = (v_{t - dt / 2} + v_{t + dt / 2}) / 2.
    p.setVelocity(tempVelocity); 

    /*if (debug) {
      cout << "//===========================================//" << endl
        << "// Particle Index: " << i << "   Num: " << i+1 << endl
        << particles[i] << endl;
    }*/
  }
}

// initialize. v_{-dt/2} = v_{0} - a_{0} * dt / 2
void ParticleSystem::initializeLeapFrog(const float& dt) {
  for(unsigned int i = 0; i < particles.size(); i++) {
    particles[i].setVelocityHalf(particles[i].getVelocity() - ((this->grav + particles[i].getAcceleration()) * dt) / 2.0f);
  }
}


void ParticleSystem::checkBoundary(Point3D* position, Vector* velocity, Vector* velocityHalf) {
  if (!position || !velocity || !velocityHalf) return;
  // temporary variables just so it will compile. these define the "boundaries" of box
  // if goes past boundaries, reflect back.
  if(position->getX() > MAX_X) {
    position->setX(MAX_X);
    bouncebackVelocity(velocity, Vector(-1, 0, 0));
    bouncebackVelocity(velocityHalf, Vector(-1, 0, 0));
	velocity->setX(-REST_COEFF * velocity->getX());
    velocityHalf->setX(-REST_COEFF * velocityHalf->getX());
  }
  else if (position->getX() < MIN_X) {
    position->setX(MIN_X);
    bouncebackVelocity(velocity, Vector(1, 0, 0));
    bouncebackVelocity(velocityHalf, Vector(1, 0, 0));
    velocity->setX(-REST_COEFF * velocity->getX());
    velocityHalf->setX(-REST_COEFF * velocityHalf->getX());
  }
  if(position->getY() > MAX_Y) {
    position->setY(MAX_Y);
    bouncebackVelocity(velocity, Vector(0, -1, 0));
    bouncebackVelocity(velocityHalf, Vector(0, -1, 0));
    velocity->setY(-REST_COEFF * velocity->getY());
    velocityHalf->setY(-REST_COEFF * velocityHalf->getY());
  }
  else if (position->getY() < MIN_Y) {
    position->setY(MIN_Y);
    bouncebackVelocity(velocity, Vector(0, 1, 0));
    bouncebackVelocity(velocityHalf, Vector(0, 1, 0));
    velocity->setY(-REST_COEFF * velocity->getY());
    velocityHalf->setY(-REST_COEFF * velocityHalf->getY());
  }
  if(position->getZ() > MAX_Z) {
    position->setZ(MAX_Z);
    bouncebackVelocity(velocity, Vector(0, 0, -1));
    bouncebackVelocity(velocityHalf, Vector(0, 0, -1));
    velocity->setZ(-REST_COEFF* velocity->getZ());
    velocityHalf->setZ(-REST_COEFF * velocityHalf->getZ());
  }
  else if (position->getZ() < MIN_Z) {
    position->setZ(MIN_Z);
    bouncebackVelocity(velocity, Vector(0, 0, 1));
    bouncebackVelocity(velocityHalf, Vector(0, 0, 1));
    velocity->setZ(-REST_COEFF* velocity->getZ());
    velocityHalf->setZ(-REST_COEFF * velocityHalf->getZ());
  }
}

void ParticleSystem::bouncebackVelocity(Vector* velocity, Vector normal) {

  //*velocity = *velocity - (normal * (*velocity).dotProduct(normal) * (1 + REST_COEFF));
}

void ParticleSystem::addParticle(Particle& p) {
  particles.push_back(p);
  grid.addParticle(p);
}
