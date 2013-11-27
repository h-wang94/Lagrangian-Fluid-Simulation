#include <iostream>
#include "Particle.h"

Particle::Particle(float mass, float volume, float pressure, float coeffVis, Vector velocity, Point3D position) {
  this->mass = mass;
  this->volume = volume;
  this->pressure = pressure;
  this->velocity = velocity;
  this->position = position;
}

float Particle::getMass() const {
  return this->mass;
}

float Particle::getVolume() const {
  return this->volume;
}

float Particle::getPressure() const {
  return this->pressure;
}

Vector Particle::getVelocity() const {
  return this->velocity;
}

Point3D Particle::getPosition() const {
  return this->position;
}

void Particle::setMass(float mass) {
  this->mass = mass;
}

void Particle::setVolume(float volume) {
  this->volume = volume;
}
void Particle::setPressure(float pressure) {
  this->pressure = pressure;
}

void Particle::setVelocity(Vector velocity) {
  this->velocity = velocity;
}

void Particle::setPosition(Point3D position) {
  this->position = position;
}

std::ostream& operator <<(ostream& outs, const Particle& particle) {
  outs << "//==================================================//" << endl
       << "// Mass: " << particle.getMass() << endl
       << "// Volume: " << particle.getVolume() << endl
       << "// Pressure: " << particle.getPressure() << endl
       << "// Velocity: " << particle.getVelocity() << endl
       << "// Position: " << particle.getPosition() << endl
       << "//==================================================//" << endl;
  return outs;
}
