#ifndef _Point3D_h_
#define _Point3D_h_
//Point3D class
//By Elliott Ison
#include "Vector.h"
#include <iostream>
#include <iomanip>
using namespace std;

class Point3D {

    /*
    Point3D consists of a position at (x, y, z)
    */

    public:
        Point3D();
        Point3D(float x, float y, float z);

        float getX() const;
        float getY() const;
        float getZ() const;
        void setX(float n);
        void setY(float n);
        void setZ(float n);
        Vector getDifference(Point3D& p);
        Vector pToVector(Point3D& p);
        Vector operator -(const Point3D& rhs) const;
        const Point3D operator +(const Point3D& rhs) const;
        const Point3D operator +(const Vector& rhs) const;
        const Point3D operator *(const float& scalar) const;

    private:
        float x;
        float y;
        float z;

};

ostream& operator <<(ostream& outs, const Point3D& point);
// end class
#endif
