/*
 * Ray.h
 *
 *  Created on: Dec 7, 2013
 *      Author: Owner
 */

#ifndef RAY_H_
#define RAY_H_

#include "Point3D.h"
#include "Vector.h"

class Ray {
public:
	Ray();
	Ray(Point3D start, Vector direction);
	Point3D getOrigin();
	Vector getDirection();
	void setOrigin(const Point3D &);
	void setDirection(const Vector &);
	Point3D getPointAtT(float t);
private:
	Point3D start;
	Vector direction;
};

#endif /* RAY_H_ */
