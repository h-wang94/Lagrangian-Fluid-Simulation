#ifndef _Particle_h_
#define _Particle_h_
#include "Vector.h"
#include "Point3D.h"

class Particle {
  public:
    Particle(void) {};
    Particle(float mass, float density, float pressure, float coeffVis, Vector velocity, Vector velocityHalf, Vector acceleration, Point3D position);
    //Particle(float density, float pressure, float coeffVis, Vector velocity, Point3D position);
    ~Particle(void) {};

    // getters and setters
    float getMass() const;
    float getVolume() const;
    float getPressure() const;
    float getDensity() const;
    Vector getVelocity() const;
    Vector getVelocityHalf() const;
    Vector getAcceleration() const;
    Point3D getPosition() const;

    void setMass(float mass);
    //void setVolume(float volume);
    void setPressure(float pressure);
    void setVelocity(Vector velocity);
    void setVelocityHalf(Vector velocityHalf);
    void setAcceleration(Vector acceleration);
    void setPosition(Point3D position);
    void setDensity(float density);

  private:
    float mass;
    //float volume;
    float density; // can be calculated from mass and volume. m/v. but we need it because it derives from m and v.
    float pressure; // is pressure given in some other way?
    Vector velocity;
    Vector velocityHalf; //velocity in halfstep for leapfrog integration. velocity = at time t, velocityHalf = at time t + 1/2
    Vector acceleration; // acceleration updated through F/rho;
    Point3D position;
    float coeffVis;

};

std::ostream& operator<<(ostream& outs, const Particle& particle);

#endif
