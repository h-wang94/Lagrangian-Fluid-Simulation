#include "Vector.h"

//begin methods
Vector::Vector() {
    this->x = 0;
    this->y = 0;
    this->z = 0;
}

Vector::Vector(float x, float y, float z) {
    this->x = x;
    this->y = y;
    this->z = z;
}

float Vector::getX() const {
    return this->x;
}

float Vector::getY() const {
    return this->y;
}

float Vector::getZ() const {
    return this->z;
}

void Vector::setX(float n){
    this->x = n;
}

void Vector::setY(float n){
    this->y = n;
}

void Vector::setZ(float n){
    this->z = n;
}

float Vector::getMagnitude(){
    //magnitude is sqrt(x^2+y^2+z^2)
    return sqrt((this->x*this->x)+(this->y*this->y)+(this->z*this->z));
}

float Vector::dotProduct(Vector u){
    return (this->getX()*u.getX())+(this->getY()*u.getY())+(this->getZ()*u.getZ());
}

Vector Vector::crossProduct(Vector u){ //(this vector) x (u vector) = (t vector)
    float x = (this->getY()*u.getZ()) - (this->getZ()*u.getY());
    float y = (this->getZ()*u.getX()) - (this->getX()*u.getZ());
    float z = (this->getX()*u.getY()) - (this->getY()*u.getX());
    return Vector(x, y, z);
}

void Vector::normalize() {
  float m = this->getMagnitude();
  this->x = this->x / m;
  this->y = this->y / m;
  this->z = this->z / m;

}

Vector& Vector::operator +=(const Vector& rhs) {
  this->x += rhs.x;
  this->y += rhs.y;
  this->z += rhs.z;
  return *this;
}

Vector& Vector::operator -=(const Vector& rhs) {
  this->x -= rhs.x;
  this->y -= rhs.y;
  this->z -= rhs.z;
  return *this;
}

const Vector Vector::operator +(const Vector& rhs) const {
  return Vector(*this) += rhs;
}

const Vector Vector::operator -(const Vector& rhs) const {
  return Vector(*this) -= rhs;
}

const Vector Vector::operator *(const float scalar) const {
  Vector temp;
  temp.x = this->x * scalar;
  temp.y = this->y * scalar;
  temp.z = this->z * scalar;
  return temp;
}

const Vector Vector::operator /(const float scalar) const {
  Vector temp;
  temp.x = this->x / scalar;
  temp.y = this->y / scalar;
  temp.z = this->z / scalar;
  return temp;
}

ostream& operator <<(ostream& outs, const Vector& vector) {
  outs << "<" << vector.getX() << ", "
    << vector.getY() << ", "
    << vector.getZ() << ">";
  return outs;
}
