#include <iostream>
#include <cstdlib>
#include "Particle.h"

#define MAX_X .09
#define MAX_Y .5
#define MAX_Z .05

// Constructor for a water particle in random locations
Particle::Particle(void) {
  this->mass = 0.02f;
  this->pressure = 0.0;
  this->stiffness = 3.0;
  this->restDensity = 998.29;
  this->density = 0;
  this->viscosity = 3.5;
  float x = MAX_X*((2*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX)))-1);
  float y = MAX_Y*((2*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX)))-1)-.5;
  float z = MAX_Z*((2*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX)))-1);
  this->position = Point3D(x, y, z);
  this->oldPosition = position;
  this->supportRadius = 0.0457f;
  this->velocity = Vector(0, 0, 0);
  this->color = Vector(0, 0, .8f);
  this->opacity = 0.8f;
  this->threshold = 7.065f;
  this->surfTension = 0.0728f;
}

Water::Water(void) {
  this->mass = 0.02f;
  this->pressure = 0.0f;
  this->stiffness = 3.0f;
  this->restDensity = 998.29f;
  this->density = 0.0f;
  this->viscosity = 3.5f;
  float x = MAX_X*((2*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX)))-1);
  float y = MAX_Y*((2*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX)))-1)-.5;
  float z = MAX_Z*((2*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX)))-1);
  this->position = Point3D(x, y, z);
  this->oldPosition = position;
  this->supportRadius = 0.0457f;
  this->velocity = Vector(0,0,0);
  this->color = Vector(0, 0, .8f);
  this->opacity = 0.5f;
  this->restCoeff = 0.9;
  this->threshold = 7.065f;
  this->surfTension = 0.0728f;
}

Mucus::Mucus(void) {
  this->mass = 0.04f;
  this->pressure = 0.0f;
  this->stiffness = 5.0f;
  this->restDensity = 1000.0f;
  this->density = 0.0f;
  this->viscosity = 36.0f;
  float x = MAX_X*((2*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX)))-1);
  float y = MAX_Y*((2*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX)))-1)-.5;
  float z = MAX_Z*((2*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX)))-1);
  this->position = Point3D(x, y, z);
  this->oldPosition = position;
  this->supportRadius = 0.0726f;
  this->velocity = Vector(0, 0, 0);
  this->color = Vector(0, .8f, 0);
  this->opacity = 0.5f;
  this->restCoeff = 0.8;
  this->threshold = 5.0f;
  this->surfTension = 6.0f;
}

Particle::Particle(float mass, float pressure, float stiffness, float restDensity, float density, float viscosity, float restCoeff, float threshold, float surfTension, Point3D position, Vector velocity, Vector color, float opacity) {
  this->mass = mass;
  this->pressure = pressure;
  this->stiffness = stiffness;
  this->restDensity = restDensity;
  this->density = density;
  this->viscosity = viscosity;
  this->position = position;
  this->oldPosition = position;
  this->velocity = velocity;
  this->color = Vector(.2,0,1);
}

int Particle::getHashID() const {
	
  return this->hashID;
}

void Particle::setHashID(const int i){
  this->hashID = i;
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

Vector Particle::getColor() const {
  return this->color;
}

float Particle::getOpacity() const {
  return this->opacity;
}

float Particle::getRestCoeff() const {
  return this->restCoeff;
}

float Particle::getSupportRadius() const {
  return this->supportRadius;
}

float Particle::getThreshold() const {
  return this->threshold;
}

float Particle::getSurfTension() const {
  return this->surfTension;
}

void Particle::setMass(const float mass) {
  this->mass = mass;
}

void Particle::setDensity(const float density) {
  this->density = density;
}

void Particle::setViscosity(const float viscosity) {
  this->density = viscosity;
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

void Particle::setColor(const Vector color) {
  this->color = color;
}

void Particle::setOpacity(const float opacity) {
  this->opacity = opacity;
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
