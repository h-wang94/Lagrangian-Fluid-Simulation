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
	CubeVertex(Point3D p); //representation for a Point3D for use with ParticleSystem
	void setColor(float c);
	float getColor();
	Particle getParticle();
	Point3D getPosition();

private:
	Point3D position;
	float colorFuncVal;
};

class Cube {
public:
	Cube();
	Cube(const vector<CubeVertex> &);
	vector<CubeVertex> getVertices();
	void setVertices(const vector<CubeVertex> &);
	void getTriangles(const float &isolevel, vector<float> &);
	void updateColors(); // debugging only (probably)

//private:
	int getCutVertices();
	int getCutEdges(int v);
	vector<int> getVertexNumsFromEdge(int edge);
	vector<CubeVertex> vertices;
};




#endif /* RECONSTRUCT_H_ */
