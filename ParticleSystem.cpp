#include <iostream>
#include "ParticleSystem.h"
#include "Particle.h"

#define PI 3.14159265

//====================================================================//
// Steps to intialize SPH system
// 1. Create fluid material
// 2. Create n particles and set position, initial velocities, and mass. Fill in missing values.
// 3. Smoothing kernels to compute compact support radius
// 4. Create spatial hashing data structure(?) and insert each particle
// 5. Create collision objects
// 6. Initialize Leapfrog integrator for all particles
//====================================================================//
// Compute density and pressure
// For each particle,
// 1. Search neighborhood using spatial hashing
// 2. Compute mass density but only iterate over particles in neighborhood
// 3. Compute pressure and use material density as rest density
//====================================================================//
// Compute internal forces
// For each particle,
// 1. Search neighborhood. Repeat because all SPH forces depend on mass-density which must be computed
//    beforea ny force density.
// 2. compute pressure force density acting on the particle
// 3. compute the viscosity force density acting on the particle
// F_internal = F_pressure + F_viscosity
//====================================================================//
// Compute external forces 
//====================================================================//
// Time integration and collision handling
//====================================================================//
// Render particles
//====================================================================//
ParticleSystem::ParticleSystem(Vector grav){
	this->grav = grav;
}

void ParticleSystem::initialize() {
  this->setDensities();
  this->computePressure(1.0f, 1.0f);
  this->computeForces();
  this->initializeLeapFrog(0.5);
}

void ParticleSystem::update(float timestep){
	this->setDensities();//for each particle, compute particle's density
  this->computePressure(1.0f, 1.0f); // now compute each particle's pressure. randomly put in numbers.
	//TODO:
	//for each particle
		//evaluate net force 
		//compute accel
		//leapfrog integration
}

void ParticleSystem::computeForces(){
  // f = rho * g - pressure forces + viscosity forces
  Vector gravity, pressure, viscosity, force;
  for(unsigned int i = 0; i < particles.size(); i++) {
    // calculate gravity forces
    gravity = gravityForce(particles[i]);
    // calculate pressure forces
    pressure = pressureForce(particles[i], i);
    // calculate viscosity forces
    viscosity = viscosityForce(particles[i], i);

    force = gravity - pressure + viscosity;
    particles[i].setAcceleration(force / particles[i].getDensity()); // why density...
  }
}

void ParticleSystem::computePressure(const float stiffness, const float restDensity) {
  for(unsigned int i = 0; i < particles.size(); i++) {
    // p = k ( (p / p0)^7 - 1)
    particles[i].setPressure(stiffness * (pow((particles[i].getDensity() / restDensity), 7.0f) - 1.0f));
  }
}

void ParticleSystem::setDensities(){
	float h = 5.0;															//CHANGE LATER (smoothing distance)
	float tol = .000001;													//CHANGE LATER (tolerance to be counted as irrelevant particle)

	for(unsigned int i = 0; i < particles.size(); i++){
		float density = 0.0;

		for(unsigned int j = 0; j < particles.size(); j++){									//INEFFICIENT for now; going to find way to only take into 
																					//account particles near particle[i]

			Vector dist = particles[i].getPosition() - particles[j].getPosition();	//need distance between two particles

			float kernel = defaultKernel(dist, h);
			if(kernel > tol){
				float density = density + kernel * particles[j].getMass();			//add on to the density for particle particles[i]
			}
		}
		particles[i].setDensity(density);											//set the particle[i]'s density to particle[i]
	}
}

Vector ParticleSystem::gravityForce(Particle& p) {
  return grav * p.getDensity(); // why density...
}

Vector ParticleSystem::pressureForce(Particle& p, unsigned const int i) {
  Vector pressure;
  float coeff;
  unsigned int j;
  // this is so stupid...........but ill think of a better way. i dont wanna do checks cause it might make a difference since this is computed every time for every particle
  for(j = 0; j < i; j++) {
    coeff = (p.getPressure() + particles[j].getPressure()) / 2.0f * particles[j].getVolume();
    pressure += pressGradientKernel(p.getPosition() - particles[j].getPosition(), 1.0f) * coeff;
  }
  for(j = j + 1; j < particles.size(); j++) {
    coeff = (p.getPressure() + particles[j].getPressure()) / 2.0f * particles[j].getVolume();
    pressure += pressGradientKernel(p.getPosition() - particles[j].getPosition(), 1.0f) * coeff;
  }
  return pressure;
}

Vector ParticleSystem::viscosityForce(Particle& p, unsigned const int i) {
  Vector viscosity;
  Vector coeff;
  unsigned int j;
  for(j = 0; j < i; j++) {
    coeff = (particles[j].getVelocity() - p.getVelocity()) * particles[j].getVolume();
    viscosity += coeff * viscLaplacianKernel(p.getPosition() - particles[i].getPosition(), 1.0f);
  }
  for(j = j + 1; j < particles.size(); j++) {
    coeff = (particles[j].getVelocity() - p.getVelocity()) * particles[j].getVolume();
    viscosity += coeff * viscLaplacianKernel(p.getPosition() - particles[i].getPosition(), 1.0f);
  }
  return viscosity;
}

// Poly6 Kernels used for everything except pressure and viscosity forces 
// Notes: http://image.diku.dk/projects/media/kelager.06.pdf (Page 16)
// Less expensive compared to gaussian one because of computation of e and no square roots
// Not sure whether I'm supposed to normalize Vector r though

float ParticleSystem::defaultKernel(Vector r, const float h) {
  float rMag = r.getMagnitude();
  if (rMag >= 0 && rMag <= h) {
    return (315.0f * pow((pow(h, 2.0f) - pow(rMag, 2.0f)),3.0) / (64.0f * PI * pow(h, 9.0f)));
  } else {
    return 0;
  }
}

// gradient and laplacian of poly6 kernels. prob not needed if we use the spiky and viscosity kernels for other calculations
/*Vector ParticleSystem::gradientKernel(Vector r, float h) {
  float rMag = r.getMagnitude();
  float coeff = pow(pow(h, 2.0f) - pow(rMag, 2.0f), 2.0f) * -945 / (32 * PI * pow(h, 9.0f));
  return r * coeff;
}

float ParticleSystem::laplacianKernel(Vector r, float h) {
  float rMag = r.getMagnitude();
  return (-945 / (32 * PI * pow(h, 9.0f))) * (pow(h, 2.0f) - pow(rMag, 2.0f)) * (3 * pow(h, 2.0f) - 7 * pow(rMag, 2.0f));
}*/

// Spiky Kernel to calculate pressure 
// Give more repulsive pressure at short distance and thus avoids clustering.
Vector ParticleSystem::pressGradientKernel(Vector r, const float h) {
  float rMag = r.getMagnitude();
  float coeff = (-45 * pow((h - rMag), 2.0f)) / (PI * pow(h, 6.0f) * rMag);
  return r * coeff;
}

// Viscosity kernel to calculate viscosity
// Gives more stable viscosity and makes it possible to damp simulation better
// Laplacian in poly6 kernel becomes negative really fast. The viscosity kernel's laplacian is positive everywhere
float ParticleSystem::viscLaplacianKernel(Vector r, const float h) {
  float rMag = r.getMagnitude();
  return (45 * (h - rMag)) / (PI * pow(h, 6.0f));

}

// Leap frog integration. Takes in a dt and will loop through all the particles in our system.
// Has not been tested.
// Can be changed to work with leapfrogging just a certain particle.
void ParticleSystem::leapFrog(const float dt) {
  for(unsigned int i = 0; i < particles.size(); i++) {
    Particle& p = particles[i];
    // get velocity at time t - dt/2. v_{t - dt/2}
    Vector old = p.getVelocityHalf(); 
    // velocity at time t + dt/2. v_{t + dt/2} = v_{t - dt / 2} + a * dt
    p.setVelocityHalf(old + (p.getAcceleration() * dt)); 
    // set position at time t. pos_{t + dt} = pos_{t} + v_{t + dt / 2} * dt
    p.setPosition(p.getPosition() + (p.getVelocityHalf() * dt)); 
    // use midpoint approximation for velocity at time t. v_{t} = (v_{t - dt / 2} + v_{t + dt / 2}) / 2.
    p.setVelocity((old + p.getVelocityHalf()) / 2.0f); 
  }
}

// initialize. v_{-dt/2} = v_{0} - a_{0} * dt / 2
void ParticleSystem::initializeLeapFrog(const float dt) {
  for(unsigned int i = 0; i < particles.size(); i++) {
    particles[i].setVelocityHalf(particles[i].getVelocity() - dt * particles[i].getAcceleration() / 2.0f);
  }
}
