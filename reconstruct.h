/*
 * reconstruct.h
 *
 *  Created on: Dec 3, 2013
 *      Author: Owner
 */

#ifndef RECONSTRUCT_H_
#define RECONSTRUCT_H_

#include "ParticleSystem.h"

extern ParticleSystem pSystem;

class CubeVertex {
public:
	CubeVertex();
	CubeVertex(Particle p); //representation for a Point3D for use with ParticleSystem
	void setColor(float c);
	float getColor();
	Particle getParticle();

private:
	Particle particle;
	float colorFuncVal;
};

class Cube {
public:
	Cube();
	Cube(const vector<CubeVertex> &);
	vector<CubeVertex> getVertices();
	void setVertices(const vector<CubeVertex> &);
	vector<Point3D> getTriangles(const float &isolevel);

private:
	__int8 getCutVertices();
	int getCutEdges(__int8 v);
	vector<int> getVertexNumsFromEdge(int edge);
	vector<CubeVertex> vertices;
};




#endif /* RECONSTRUCT_H_ */
