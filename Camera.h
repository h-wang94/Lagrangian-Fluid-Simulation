/*
 * Camera.h
 *
 *  Created on: Dec 7, 2013
 *      Author: Owner
 */

#ifndef CAMERA_H_
#define CAMERA_H_

#include "Ray.h"

class Camera {
	friend class Sampler;
public:
	Camera();
	Camera(float imgH, float imgW, float fovV, Point3D lookFrom, Point3D lookAt, Vector up);
private:
	float imgHeight, imgWidth, R, T;
	Vector xAxis, yAxis, zAxis;
	Point3D lookFrom;
};


#endif /* CAMERA_H_ */
