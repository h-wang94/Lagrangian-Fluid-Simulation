#include <iostream>
#include "ParticleSystem.h"
#include "Particle.h"

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