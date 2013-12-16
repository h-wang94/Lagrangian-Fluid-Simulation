#include "ParticleSystem.h"
#include <omp.h>
#include <time.h>

#define PI 3.14159265
#define MAX_NAIVE 50

ParticleSystem::ParticleSystem() {
  this->grav = Vector(0,0,-9.8);
  setBoundaries(-.1, -1, -.1, .1, 1, .1);
}

ParticleSystem::ParticleSystem(Vector grav){
  this->grav = grav;
  this->numRowBoxes = 100;
  this->grid = SpatialGrid(0.0457);
  setBoundaries(-.1, -1, -.1, .1, 1, .1);
}

ParticleSystem::ParticleSystem(Vector grav, float min_x, float min_y, float min_z, float max_x, float max_y, float max_z) {
  this->grav = grav;
  this->numRowBoxes = 100;
  this->grid = SpatialGrid(0.0457);
  setBoundaries(min_x, min_y, min_z, max_x, max_y, max_z);
}

void ParticleSystem::setBoundaries(float min_x, float min_y, float min_z, float max_x, float max_y, float max_z) {
  this->MIN_X = min_x;
  this->MIN_Y = min_y;
  this->MIN_Z = min_z;
  this->MAX_X = max_x;
  this->MAX_Y = max_y;
  this->MAX_Z = max_z;
}

void ParticleSystem::initialize(float timestep) {
  this->setDensities();
  this->computePressure();
  this->computeForces();
  this->initializeLeapFrog(-timestep);
}

void ParticleSystem::update(float timestep){
  this->setDensities();//for each particle, compute particle's density
  this->computePressure(); // now compute each particle's pressure
  this->computeForces(); // compute forces and set acceleration of the particles
  this->leapFrog(timestep); // utilize leapfrog integration to figure out position and new velocity
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
  Vector gravity, force, total, tension;
#pragma omp parallel for firstprivate(gravity, total, force, tension)
  for(unsigned int i = 0; i < particles.size(); i++) {
    gravity = gravityForce(particles[i]);
    total = press_visc(particles[i], i); // viscosity - pressure forces

    force = gravity + total;
    Vector acceleration = force / particles[i].getDensity();
    particles[i].setAcceleration(acceleration);
  }
}

void ParticleSystem::computePressure() {
  float pressure = 0;
  unsigned int i = 0;
#pragma omp parallel for firstprivate(pressure)
  // loop unroll by 4
  for(i = 0; i < particles.size()/4*4; i+=4) {
    // Equation: p = k ( (p / p0)^7 - 1)
    pressure = particles[i].getStiffness() * (pow((particles[i].getDensity() / particles[i].getRestDensity()), 7.0f) - 1.0f);
    particles[i].setPressure(pressure);
    pressure = particles[i+1].getStiffness() * (pow((particles[i+1].getDensity() / particles[i+1].getRestDensity()), 7.0f) - 1.0f);
    particles[i+1].setPressure(pressure);
    pressure = particles[i+2].getStiffness() * (pow((particles[i+2].getDensity() / particles[i+2].getRestDensity()), 7.0f) - 1.0f);
    particles[i+2].setPressure(pressure);
    pressure = particles[i+3].getStiffness() * (pow((particles[i+3].getDensity() / particles[i+3].getRestDensity()), 7.0f) - 1.0f);
    particles[i+3].setPressure(pressure);
  }
  for(; i < particles.size(); i++) {
    pressure = particles[i].getStiffness() * (pow((particles[i].getDensity() / particles[i].getRestDensity()), 7.0f) - 1.0f);
    particles[i].setPressure(pressure);
  }
}

void ParticleSystem::setDensities(){
  float density = 0;
  std::vector<Particle> list;

  std::map<int, std::vector<int> > hashMap;

  for(unsigned int i = 0; i < particles.size(); i++){//create the hash map (O(n) time)

	Point3D pos = particles[i].getPosition();
	float x = pos.getX();
	float y = pos.getY();
	float z = pos.getZ();

	int xindex = (int)floor((x+10)/(.0457));
	int yindex = (int)floor((y+10)/(.0457));
	int zindex = (int)floor((z+10)/(.0457));

	int hash = (xindex * 10619863)^(yindex * 94418953)^(zindex * 54018521);

	hashMap[hash].push_back(i);

  }

//#pragma omp parallel for firstprivate(density, list)
  for(unsigned int i = 0; i < particles.size(); i++){
    density = 0;

    if(particles.size() > MAX_NAIVE){//doing spatial hashing
      //list = grid.getNeighbors(particles[i]);

		Point3D pos = particles[i].getPosition();
		float x = pos.getX();
		float y = pos.getY();
		float z = pos.getZ();

		int xindex = (int)floor((x+10)/(.0457));
		int yindex = (int)floor((y+10)/(.0457));
		int zindex = (int)floor((z+10)/(.0457));

		int hash = (xindex * 10619863)^(yindex * 94418953)^(zindex * 54018521);
		std::vector<int> ints = hashMap[hash];
#pragma omp parallel for firstprivate(c)
		for(int c = 0; c < ints.size(); c++){//get surrounding parts
			list.push_back(particles[ints[c]]);
		}
		for(int a = -1; a <= 1; a++){
			for(int b = -1; b <= 1; b++){
				for(int c = -1; c <= 1; c++){//get surrounding boxes 
					int hash2 = ((xindex+a) * 10619863)^((yindex+b) * 94418953)^((zindex+c) * 54018521);
					std::vector<int> ints2 = hashMap[hash2];
					#pragma omp parallel for firstprivate(d)
					for(int d = 0; d < ints2.size(); d++){
						list.push_back(particles[ints2[d]]);
					}
				}
			}
		}
    for(unsigned int j = 0; j < list.size(); j++){ //loop through neighbors
        Vector diff = particles[i].getPosition() - list[j].getPosition();
        if (diff.getMagnitude() <= particles[i].getSupportRadius()) {
        density += defaultKernel(diff, particles[i].getSupportRadius()) * list[j].getMass();
        }
      }
	list.clear();
    }

    else{
      for(unsigned int j = 0; j < particles.size(); j++){ 
        Vector diff = particles[i].getPosition() - particles[j].getPosition();
        if (diff.getMagnitude() <= particles[i].getSupportRadius()) {
        density += defaultKernel(diff, particles[i].getSupportRadius()) * particles[j].getMass();
        }
      }
    }
    if(density == 0){
      density = particles[i].getMass() / .00000001;
    }
    particles[i].setDensity(density);
  }
}

Vector ParticleSystem::gravityForce(Particle& p) {
  return grav * p.getDensity(); 
}

Vector ParticleSystem::press_visc(Particle& p , unsigned const int& i) {
  unsigned int j = 0;
  Vector pressure, viscosity, diff;
  float pCoeff;
  Vector vCoeff;
  for(j = 0; j < i; j++) {
    diff = p.getPosition() - particles[j].getPosition();
    if(diff.getMagnitude() <= p.getSupportRadius()) {
      pCoeff = (p.getPressure() + particles[j].getPressure()) / 2.0 * particles[j].getVolume();
      pressure += pressGradientKernel(diff, p.getSupportRadius()) * pCoeff;
      vCoeff = (particles[j].getVelocity() - p.getVelocity()) * particles[j].getVolume();
      viscosity += vCoeff * viscLaplacianKernel(diff, p.getSupportRadius());
    }
  }
  for(j = j + 1; j < particles.size(); j++) {
    diff = p.getPosition() - particles[j].getPosition();
    if(diff.getMagnitude() <= p.getSupportRadius()) {
      pCoeff = (p.getPressure() + particles[j].getPressure()) / 2.0 * particles[j].getVolume();
      pressure += pressGradientKernel(diff, p.getSupportRadius()) * pCoeff;
      vCoeff = (particles[j].getVelocity() - p.getVelocity()) * particles[j].getVolume();
      viscosity += vCoeff * viscLaplacianKernel(diff, p.getSupportRadius());
    }
  }
  viscosity = viscosity * p.getViscosity();
  return viscosity - pressure;
}

Vector ParticleSystem::pressureForce(Particle& p, unsigned const int& i) {
  Vector pressure;
  float coeff = 0;
  unsigned int j = 0;
  // this is so stupid...........but ill think of a better way. i dont wanna do checks cause it might make a difference since this is computed every time for every particle
  for(j = 0; j < i; j++) {
    Vector diff = p.getPosition() - particles[j].getPosition();
    if (diff.getMagnitude() <= p.getSupportRadius()) {
      coeff = (p.getPressure() + particles[j].getPressure()) / 2.0f * particles[j].getVolume();
      pressure += pressGradientKernel(diff, p.getSupportRadius()) * coeff;
    }
  }
  for(j = j + 1; j < particles.size(); j++) {
    Vector diff = p.getPosition() - particles[j].getPosition();
    if (diff.getMagnitude() <= p.getSupportRadius()) {
      coeff = (p.getPressure() + particles[j].getPressure()) / 2.0f * particles[j].getVolume();
      pressure += pressGradientKernel(diff, p.getSupportRadius()) * coeff;
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
    if (diff.getMagnitude() <= p.getSupportRadius()) {
      coeff = (particles[j].getVelocity() - p.getVelocity()) * particles[j].getVolume();
      viscosity += coeff * viscLaplacianKernel(diff, p.getSupportRadius());
    }
  }
  for(j = j + 1; j < particles.size(); j++) {
    Vector diff = p.getPosition() - particles[j].getPosition();
    if (diff.getMagnitude() <= p.getSupportRadius()) {
      coeff = (particles[j].getVelocity() - p.getVelocity()) * particles[j].getVolume();
      viscosity += coeff * viscLaplacianKernel(diff, p.getSupportRadius());
    }
  }
  viscosity = viscosity * p.getViscosity(); 
  return viscosity;
}

Vector ParticleSystem::tensionForce(Particle& p, unsigned const int& i) {
  Vector normal = surfaceNormal(p, i);
  if (normal.getMagnitude() > p.getThreshold()) { // threshold
    normal.normalize();
    return normal * curvature(p, i) * (p.getSurfTension()); // surface tension constant for water
  }
  return Vector(0, 0, 0);
}

float ParticleSystem::colorFunction(Particle& p) {
  float color = 0;
  Vector diff;
  unsigned int j;
  for(j = 0; j < particles.size(); j++) {
    diff = p.getPosition() - particles[j].getPosition();
    if (diff.getMagnitude() <= p.getSupportRadius()) {
      color += particles[j].getMass() / particles[j].getDensity() * defaultKernel(diff, p.getSupportRadius());
    }
  }
  return color;
}

vector<Particle> ParticleSystem::getNeighbors(Particle &p) {
  Vector diff;
  vector<Particle> neighbors;
  unsigned int j;
  for(j = 0; j < particles.size(); j++) {
    diff = p.getPosition() - particles[j].getPosition();
    if (diff.getMagnitude() <= p.getSupportRadius()) {
      neighbors.push_back(particles[j]);
    }
  }
  return neighbors;
}

Vector ParticleSystem::surfaceNormal(Particle& p, unsigned const int& i) {
  Vector normal = Vector(0, 0, 0);
  Vector diff;
  unsigned int j;
  for(j = 0; j < i; j++) {
    diff = p.getPosition() - particles[j].getPosition();
    if (diff.getMagnitude() <= p.getSupportRadius()) {
      normal += gradientKernel(diff, p.getSupportRadius()) * particles[j].getMass() / particles[j].getDensity();
    }
  }
  for(j = j + 1; j < particles.size(); j++) {
    diff = p.getPosition() - particles[j].getPosition();
    if (diff.getMagnitude() <= p.getSupportRadius()) {
      normal += gradientKernel(diff, p.getSupportRadius()) * particles[j].getMass() / particles[j].getDensity();
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
    if (diff.getMagnitude() <= p.getSupportRadius()) {
      curvature += particles[j].getMass() / particles[j].getDensity() * laplacianKernel(diff, p.getSupportRadius());
    }
  }
  for(j = j + 1; j < particles.size(); j++) {
    diff = p.getPosition() - particles[j].getPosition();
    if (diff.getMagnitude() <= p.getSupportRadius()) {
      curvature += particles[j].getMass() / particles[j].getDensity() * laplacianKernel(diff, p.getSupportRadius());
    }
  }
  return curvature;
}


// Smoothing distance h is half of the difference between particle i's most distant nearest neighbor and i
// Poly6 Kernels used for everything except pressure and viscosity forces 
// Notes: http://image.diku.dk/projects/media/kelager.06.pdf (Page 16)
// Less expensive compared to gaussian one because of computation of e and no square roots
// Not sure whether I'm supposed to normalize Vector r though

float ParticleSystem::defaultKernel(Vector r, const float h) {
  float rMag = r.getMagnitude();
  return (315.0f * pow((pow(h, 2.0f) - pow(rMag, 2.0f)),3.0) / (64.0f * PI * pow(h, 9.0f)));
}

/* gradient and laplacian of poly6 kernels. prob not needed if we use the spiky and viscosity kernels for other calculations*/
Vector ParticleSystem::gradientKernel(Vector r, const float h) {
  float rMag = r.getMagnitude();
  float coeff = pow(pow(h, 2.0f) - pow(rMag, 2.0f), 2.0f) * -945 / (32 * PI * pow(h, 9.0f));
  return r * coeff;
}

float ParticleSystem::laplacianKernel(Vector r, const float h) {
  float rMag = r.getMagnitude();
  return (-945 / (32 * PI * pow(h, 9.0f))) * (pow(h, 2.0f) - pow(rMag, 2.0f)) * (3 * pow(h, 2.0f) - 7 * pow(rMag, 2.0f));
}

// Spiky Kernel to calculate pressure 
// Give more repulsive pressure at float distance and thus avoids clustering.
Vector ParticleSystem::pressGradientKernel(Vector r, const float h) {
  float rMag = r.getMagnitude();
  float coeff = (-45 * pow((h - rMag), 2.0f)) / (PI * pow(h, 6.0f));
  //r.normalize();
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
// Can be changed to work with leapfrogging just a certain particle.
void ParticleSystem::leapFrog(const float& dt) {
#pragma omp parallel for
  for(unsigned int i = 0; i < particles.size(); i++) {
    Particle& p = particles[i];
    // get velocity at time t - dt/2. v_{t - dt/2}
    Vector old = p.getVelocityHalf(); 
    // velocity at time t + dt/2. v_{t + dt/2} = v_{t - dt / 2} + a * dt
    Vector tempVelocityHalf = old + p.getAcceleration() * dt;
    // set position at time t. pos_{t + dt} = pos_{t} + v_{t + dt / 2} * dt
    Point3D tempPosition = p.getPosition() + (tempVelocityHalf * dt);
    Vector tempVelocity = (old + p.getVelocityHalf()) / 2.0f;
    checkBoundary(p, &tempPosition, &tempVelocity, &tempVelocityHalf);
    p.setOldPosition(p.getPosition());
    p.setVelocityHalf(tempVelocityHalf); 
    p.setPosition(tempPosition); 
    // use midpoint approximation for velocity at time t. v_{t} = (v_{t - dt / 2} + v_{t + dt / 2}) / 2.
    p.setVelocity(tempVelocity); 

    if(particles.size() > MAX_NAIVE){
      //p = this->grid.updateBoxes(p); // update the boxes of the particles in spatialgrid
    }
  }
}

// initialize. v_{-dt/2} = v_{0} - a_{0} * dt / 2
void ParticleSystem::initializeLeapFrog(const float& dt) {
  for(unsigned int i = 0; i < particles.size(); i++) {
    particles[i].setVelocityHalf(particles[i].getVelocity() - ((this->grav + particles[i].getAcceleration()) * dt) / 2.0f);
  }
}


void ParticleSystem::checkBoundary(Particle& p, Point3D* position, Vector* velocity, Vector* velocityHalf) {
  if (!position || !velocity || !velocityHalf) return;
  // if goes past boundaries, reflect back.
  if(position->getX() > MAX_X) {
    position->setX(MAX_X);
	
    velocity->setX(-p.getRestCoeff() * velocity->getX());
    velocityHalf->setX(-p.getRestCoeff() * velocityHalf->getX());
  }
  else if (position->getX() < MIN_X) {
    position->setX(MIN_X);

    velocity->setX(-p.getRestCoeff() * velocity->getX());
    velocityHalf->setX(-p.getRestCoeff() * velocityHalf->getX());
  }
  if(position->getY() > MAX_Y) {
    position->setY(MAX_Y);
	
    velocity->setY(-p.getRestCoeff() * velocity->getY());
    velocityHalf->setY(-p.getRestCoeff() * velocityHalf->getY());
  }
  else if (position->getY() < MIN_Y) {
    position->setY(MIN_Y);

    velocity->setY(-p.getRestCoeff() * velocity->getY());
    velocityHalf->setY(-p.getRestCoeff() * velocityHalf->getY());
  }
  if(position->getZ() > MAX_Z) {
    position->setZ(MAX_Z);

    velocity->setZ(-p.getRestCoeff()* velocity->getZ());
    velocityHalf->setZ(-p.getRestCoeff() * velocityHalf->getZ());
  }
  else if (position->getZ() < MIN_Z) {
    position->setZ(MIN_Z);

    velocity->setZ(-p.getRestCoeff()* velocity->getZ());
    velocityHalf->setZ(-p.getRestCoeff() * velocityHalf->getZ());
  }
}

void ParticleSystem::addParticle(Particle& p) {
  //long id = grid.addParticle(p);
  //p.setHashID(id);

  particles.push_back(p);
}
