/*
 * Ray.cpp
 *
 *  Created on: Dec 7, 2013
 *      Author: Owner
 */

#include "Ray.h"

Ray::Ray() {
	start = Point3D();
	direction = Vector();
}

Ray::Ray(Point3D p, Vector d) {
	start = p;
	direction = d;
}

Point3D Ray::getOrigin() {
	return start;
}

Vector Ray::getDirection() {
	return direction;
}

Point3D Ray::getPointAtT(float t) {
	return start + direction * t;
}

void Ray::setOrigin(const Point3D &p) {
	start = p;
}

void Ray::setDirection(const Vector &v) {
	direction = v;
}
