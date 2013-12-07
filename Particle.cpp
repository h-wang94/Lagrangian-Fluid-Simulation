#include <iostream>
#include <cstdlib>
#include "Particle.h"

// Constructor for a water particle in random locations
Particle::Particle(void) {
  this->mass = 0.02f;
  this->pressure = 0.0;
  this->stiffness = 3.0;
  this->restDensity = 998.29;
  this->density = 0;
  this->viscosity = 3.5;
  float x = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
  float y = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
  float z = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
  this->position = Point3D(x, y, z);
  this->oldPosition = position;
  //this->velocity = Vector(rand() % 3, rand() % 3, rand() % 3);
  this->velocity = Vector(0, 0, 0);

}

Particle::Particle(float mass, float pressure, float stiffness, float restDensity, float density, float viscosity, Point3D position, Vector velocity) {
  this->mass = mass;
  this->pressure = pressure;
  this->stiffness = stiffness;
  this->restDensity = restDensity;
  this->density = density;
  this->viscosity = viscosity;
  this->position = position;
  this->oldPosition = position;
  this->velocity = velocity;
}

float Particle::getMass() const {
  return this->mass;
}

float Particle::getVolume() const {
  return (this->mass)/(this->density);
}

float Particle::getPressure() const {
  return this->pressure;
}

float Particle::getStiffness() const {
  return this->stiffness;
}

float Particle::getRestDensity() const {
  return this->restDensity;
}

float Particle::getDensity() const {
  return this->density;
}

float Particle::getViscosity() const {
  return this->viscosity;
}

Point3D Particle::getPosition() const {
  return this->position;
}

Point3D Particle::getOldPosition() const {
  return this->oldPosition;
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

void Particle::setMass(const float mass) {
  this->mass = mass;
}

void Particle::setDensity(const float density) {
  this->density = density;
}

void Particle::setPressure(const float pressure) {
  this->pressure = pressure;
}

void Particle::setPosition(Point3D position) {
  this->position = position;
}

void Particle::setOldPosition(Point3D oldPosition) {
  this->oldPosition = oldPosition;
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
	   << "// Old Position: " << particle.getOldPosition() << endl
       << "//===========================================//" << endl;
  return outs;
}
