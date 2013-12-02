#include <iostream>
#include "ParticleSystem.h"
#include "Particle.h"

#define PI 3.14159265

ParticleSystem::ParticleSystem(float grav){
	this->grav = grav;
}

void ParticleSystem::computeForces(){

}

void ParticleSystem::update(float timestep){
	this->setDensities();//for each particle, compute particle's density

	//TODO:
	//for each particle
		//evaluate net force 
		//compute accel
		//leapfrog integration
}

void ParticleSystem::setDensities(){
	
	for(int i = 0; i < particles.size(); i++){
		float density = 0.0;
		for(int j = 0; j < particles.size(); j++){									//INEFFICIENT for now; going to find way to only take into 
																					//account particles near particle[i]

			Vector dist = particles[i].getPosition() - particles[j].getPosition();	//need distance between two particles
			float d = dist.getMagnitude();											//get the distance
			float h = 5.0;															//CHANGE LATER (smoothing distance)
			float tol = .000001;													//CHANGE LATER (tolerance to be counted as irrelevant particle)
			float kernal = (1/(pow(3.14, (3/2))*pow(h,3)))*exp((d*d)/(h*h));		//Gaussian kernal smoothing function
			if(kernal > tol){
				float density = density + kernal * particles[j].getMass();			//add on to the density for particle particles[i]
			}


		}
		particles[i].setDensity(density);											//set the particle[i]'s density to particle[i]
	}
}

// Poly6 Kernels used for everything except pressure and viscosity forces 
// Notes: http://image.diku.dk/projects/media/kelager.06.pdf (Page 16)
// Less expensive compared to gaussian one because of computation of e and no square roots
// Not sure whether I'm supposed to normalize Vector r though

float ParticleSystem::defaultKernel(Vector r, float h) {
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
Vector ParticleSystem::pressGradientKernel(Vector r, float h) {
  float rMag = r.getMagnitude();
  float coeff = (-45 * pow((h - rMag), 2.0f)) / (PI * pow(h, 6.0f) * rMag);
  return r * coeff;
}

// Viscosity kernel to calculate viscosity
// Gives more stable viscosity and makes it possible to damp simulation better
// Laplacian in poly6 kernel becomes negative really fast. The viscosity kernel's laplacian is positive everywhere
float ParticleSystem::viscLaplacianKernel(Vector r, float h) {
  float rMag = r.getMagnitude();
  return (45 * (h - rMag)) / (PI * pow(h, 6.0f));

}
