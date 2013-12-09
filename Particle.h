#ifndef _Particle_h_
#define _Particle_h_
#include "Vector.h"
#include "Point3D.h"

class Particle {
  public:
    Particle(void);
    Particle(float mass, float pressure, float stiffness, float restDensity, float density, float viscosity, Point3D position, Vector velocity, Vector color=Vector(0.0f, 0.0f, 0.8f), float opacity=0.5f);
    ~Particle(void) {};

    // getters and setters
    float getMass() const;
    float getVolume() const;
    float getPressure() const;
    float getStiffness() const;
    float getRestDensity() const;
    float getDensity() const;
    float getViscosity() const;
    Point3D getPosition() const;
    Point3D getOldPosition() const;
    Vector getVelocity() const;
    Vector getVelocityHalf() const;
    Vector getAcceleration() const;
    Vector getColor() const;
    float getOpacity() const;

    void setMass(const float mass);
    void setDensity(const float density);
    void setPressure(const float pressure);
    void setPosition(const Point3D position);
    void setOldPosition(const Point3D oldPosition);
    void setVelocity(const Vector velocity);
    void setVelocityHalf(const Vector velocityHalf);
    void setAcceleration(const Vector acceleration);
    void setColor(const Vector color);
    void setOpacity(const float opacity);

  private:
    float mass;
    float pressure;
    float stiffness;
    float restDensity;
    float density; // can be calculated from mass and volume. m/v. but we need it because it derives from m and v.
    float viscosity;
    Point3D position;
    Point3D oldPosition;
    Vector velocity;
    Vector velocityHalf; //velocity in halfstep for leapfrog integration. velocity = at time t, velocityHalf = at time t - 1/2
    Vector acceleration; // acceleration updated through F/rho;
    float coeffVis; // what is this for?
    Vector color;
    float opacity;

};

std::ostream& operator<<(ostream& outs, const Particle& particle);

#endif
