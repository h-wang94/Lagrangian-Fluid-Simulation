/*
 * mesh.h
 *
 *  Created on: Dec 11, 2013
 *      Author: Owner
 */

#ifndef MESH_H_
#define MESH_H_

#include "reconstruct.h"

extern int frameWidth, frameHeight;

class Mesh {
public:
	Mesh();
	Mesh(int width, int height, int depth, float leftbound, float rightbound, float highbound, float lowbound, float closebound, float farbound);
	CubeVertex getVertexAt(int x, int y, int z);
	void setVertexAt(int x, int y, int z, const CubeVertex &c);
	void updateColors();
	void marchingCubes(vector<float> &triangles, vector<float> &normals); //returns vectors of vectors of points for triangles

private:
	std::vector<CubeVertex> mesh;
	int width, height, depth;
	float threshold;
};

#endif /* MESH_H_ */
