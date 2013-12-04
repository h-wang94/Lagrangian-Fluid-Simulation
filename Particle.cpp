#include <iostream>
#include <cstdlib>
#include "Particle.h"

Particle::Particle(void) {
  this->mass = 0.02f;
  this->pressure = 3.0;
  this->viscosity = 3.5;
  this->velocity = Vector(rand() % 10, rand() % 10, rand() % 10);
  this->position = Point3D(rand() % 10, rand() % 10, rand() % 10);

}

Particle::Particle(float mass, float pressure, float viscosity, float coeffVis, Vector velocity, Point3D position) {
  this->mass = mass;
  //this->volume = volume; //a volume refers to a group of particles not a single particle
  this->pressure = pressure;
  this->viscosity = viscosity;
  this->coeffVis = coeffVis; // not sure what this is
  this->velocity = velocity;
  //this->velocityHalf = Vector(0,0,0); // need acceleration to compute the initial v_half
  //this->acceleration = Vector(0,0,0);
  this->position = position;
  //this->density = density;//density is dependent on the amount of particles relative to this particle
}


float Particle::getMass() const {
  return this->mass;
}

float Particle::getVolume() const {
  return (this->mass)/(this->density);
}

float Particle::getViscosity() const {
  return this->viscosity;
}

float Particle::getDensity() const {
  return this->density;
}

float Particle::getPressure() const {
  return this->pressure;
}

Vector Particle::getVelocity() const {
  return this->velocity;
}

Vector Particle::getVelocityHalf() const {
  return this->velocityHalf;
}

Vector Particle::getAcceleration() const {
  return this->acceleration;
}

Point3D Particle::getPosition() const {
  return this->position;
}

void Particle::setMass(float mass) {
  this->mass = mass;
}

/*void Particle::setVolume(float volume) {
  this->volume = volume;
}*/

void Particle::setViscosity(float viscosity) {
  this->viscosity = viscosity;
}

void Particle::setDensity(float density) {
  this->density = density;
}

void Particle::setPressure(float pressure) {
  this->pressure = pressure;
}

void Particle::setVelocity(Vector velocity) {
  this->velocity = velocity;
}

void Particle::setVelocityHalf(Vector velocityHalf) {
  this->velocityHalf = velocityHalf;
}

void Particle::setAcceleration(Vector acceleration) {
  this->acceleration = acceleration;
}

void Particle::setPosition(Point3D position) {
  this->position = position;
}

std::ostream& operator <<(ostream& outs, const Particle& particle) {
  outs << "//===========================================//" << endl
       << "// Mass: " << particle.getMass() << endl
       << "// Volume: " << particle.getVolume() << endl
       << "// Density: " << particle.getDensity() << endl
       << "// Pressure: " << particle.getPressure() << endl
       << "// Velocity: " << particle.getVelocity() << endl
       << "// Velocity Half: " << particle.getVelocityHalf() << endl
       << "// Acceleration: " << particle.getAcceleration() << endl
       << "// Position: " << particle.getPosition() << endl
       << "//===========================================//" << endl;
  return outs;
}
