/*
 * Camera.cpp
 *
 *  Created on: Dec 7, 2013
 *      Author: Owner
 */

#include "Camera.h"

Camera::Camera() {
	lookFrom = Point3D(0, 0, 0);
	xAxis = Vector(1, 0, 0);
	yAxis = Vector(0, 1, 0);
	zAxis = Vector(0, 0, 1);
	imgWidth = 640;
	imgHeight = 480;
	T = .5;
	R = T * imgWidth/imgHeight;
}

Camera::Camera(float imgW, float imgH, float fovV, Point3D lookFrom, Point3D lookAt, Vector up) {
	this->lookFrom = lookFrom;
	Vector temp = (lookAt - lookFrom);
	temp.normalize();
	zAxis = temp;
	xAxis = zAxis.crossProduct(up);
	xAxis.normalize();
	yAxis = xAxis.crossProduct(zAxis);
	yAxis.normalize();
	imgHeight = imgH;
	imgWidth = imgW;
	T = tan(fovV * .5);
	R = T * imgW/imgH;
}
