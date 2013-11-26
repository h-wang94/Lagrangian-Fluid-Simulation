#include "math.h"
#ifndef _Vector_h_
#define _Vector_h_
#include <iostream>
using namespace std;
class Vector {

    /*
    Vector consists of a direction (x, y, z) (this is not a 
    unit vector!) and can return a calculated magnitude
    */

    public:
        Vector();
        Vector(float x, float y, float z);
        float getX() const;
        float getY() const;
        float getZ() const;
        void setX(float n);
        void setY(float n);
        void setZ(float n);
        float getMagnitude();
        float dotProduct(Vector u);
        Vector crossProduct(Vector u);
        void normalize();

        // operator overloading
        Vector& operator +=(const Vector& rhs);
        Vector& operator -=(const Vector& rhs);
        const Vector operator +(const Vector& rhs) const;
        const Vector operator -(const Vector& rhs) const; 
        const Vector operator *(const float scalar) const;
        const Vector operator /(const float scalar) const;

    private:
        float x;
        float y;
        float z;
};

ostream& operator <<(ostream& outs, const Vector& vector);

#endif
