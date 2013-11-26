#include "Point3D.h"

//begin methods
Point3D::Point3D(){

  this->x = 0;
  this->y = 0;
  this->z = 0;
}

Point3D::Point3D(float x, float y, float z){

  this->x = x;
  this->y = y;
  this->z = z;
}

float Point3D::getX() const {
  return this->x;
}

float Point3D::getY() const {
  return this->y;
}

float Point3D::getZ() const {
  return this->z;
}

void Point3D::setX(float n){
  this->x = n;
}

void Point3D::setY(float n){
  this->y = n;
}

void Point3D::setZ(float n){
  this->z = n;
}

Vector Point3D::getDifference(Point3D& p){//this is the starting Point3D, p is the ending Point3D.
  return Vector(p.x - this->x, p.y - this->y, p.z - this->z);
}

Vector Point3D::pToVector(Point3D& point) {
  return Vector(point.getX(), point.getY(), point.getZ());
}

Vector Point3D::operator -(const Point3D& rhs) const {
  return Vector(this->x - rhs.x, this->y - rhs.y, this->z - rhs.z);
}

const Point3D Point3D::operator +(const Point3D& rhs) const {
  return Point3D(this->x + rhs.x, this->y + rhs.y, this->z + rhs.z);
}

const Point3D Point3D::operator *(const float& scalar) const {
  return Point3D(this->x * scalar, this->y * scalar, this->z * scalar);
}

ostream& operator <<(ostream& outs, const Point3D& point) {
  outs << "(" << point.getX() << ", "
    << point.getY() << ", "
    << point.getZ() << ")";
  return outs;
}
