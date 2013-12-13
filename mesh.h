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
	Mesh(int width, int height, int depth);
	CubeVertex getVertexAt(int x, int y, int z);
	void setVertexAt(int x, int y, int z, const CubeVertex &c);
	void updateColors();
	std::vector<std::vector<Point3D> > marchingCubes(); //returns vectors of vectors of points for triangles


	std::vector<CubeVertex> mesh;
	int width, height, depth;
};

#endif /* MESH_H_ */
