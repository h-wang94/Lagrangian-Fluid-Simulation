/*
 * reconstruct.h
 *
 *  Created on: Dec 3, 2013
 *      Author: Owner
 */

#ifndef RECONSTRUCT_H_
#define RECONSTRUCT_H_

#include "ParticleSystem.h"

class Cube {
public:
	Cube();
	Cube(const vector<Point3D> &);
	vector<Point3D> getVertices();
	void setVertices(const vector<Point3D> &);
	vector<Point3D> getTriangles(const double &isolevel, const __int8 &vertices);
	float colorFunction(const ParticleSystem &ps, const vector<Particle> &neighbors, const float neighborRadius, const Point3D &position);

private:
	__int8 getCutVertices();
	int getCutEdges(__int8 v);
	int* getTriangleVertexList(__int8 v);
	int* getVertexNumsFromEdge(int edge);
	vector<Point3D> vertices;
};




#endif /* RECONSTRUCT_H_ */
