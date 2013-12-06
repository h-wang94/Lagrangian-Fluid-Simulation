/*
 * reconstruct.h
 *
 *  Created on: Dec 3, 2013
 *      Author: Owner
 */

#ifndef RECONSTRUCT_H_
#define RECONSTRUCT_H_

#include "ParticleSystem.h"

class Surface {
public:
	float colorFunction(const ParticleSystem &ps, const vector<Particle> &neighbors, const float neighborRadius, const Point3D &position);
	//float signedDistField(float position);
	//float weightAtParticle(Particle p);
};

class Cube {
public:
	Cube();
	Cube(const vector<Point3D> &);
	vector<Point3D> getVertices();
	void setVertices(const vector<Point3D> &);

	__int8 getCutVertices(Surface s);
	int getCutEdges(__int8 v);

private:
	vector<Point3D> vertices;
};




#endif /* RECONSTRUCT_H_ */
