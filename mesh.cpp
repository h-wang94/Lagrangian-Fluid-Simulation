/*
 * mesh.cpp
 *
 *  Created on: Dec 11, 2013
 *      Author: Owner
 */

#include "mesh.h"
#include <math.h>

Mesh::Mesh() {
	width = 0;
	height = 0;
	depth = 0;
	threshold = 0;
}

/* A mesh WIDTH boxes wide, HEIGHT boxes tall, DEPTH boxes deep. */
Mesh::Mesh(int width, int height, int depth, float rightbound, float leftbound, float highbound, float lowbound, float closebound, float farbound) {
	this->width = width;
	this->height = height;
	this->depth = depth;
	float cubeW = (rightbound - leftbound) / (float) width;
	float cubeH = (highbound - lowbound) / (float) height; //height of a cube
	float cubeD = (closebound - farbound) / (float) depth;

	threshold = .1 * min(cubeW, min(cubeH, cubeD));
	for (float z = closebound; z >= (farbound - threshold); z-=cubeD) {
		for (float y = highbound; y >= (lowbound - threshold); y-=cubeH) {
			for (float x = leftbound; x <= (rightbound + threshold); x+=cubeW) {
				mesh.push_back(CubeVertex(Point3D(x, y, z)));
			}
		}
	}
}

CubeVertex Mesh::getVertexAt(int x, int y, int z) {
	return mesh[(z * (width+1) * (height+1)) + (y * (width+1)) + x]; //wat
}

void Mesh::setVertexAt(int x, int y, int z, const CubeVertex &c) {
	mesh[(z * (width+1) * (height+1)) + (y * (height+1)) + x] = c;
}

void Mesh::updateColors() {
    #pragma omp parallel for
	for (int i = 0; i < (int) mesh.size(); i++) {
		Particle p = mesh[i].getParticle();
		mesh[i].setColor(pSystem.colorFunction(p));
	}
}

void Mesh::marchingCubes(vector<float> &triangles, vector<float> &normals) {
	Cube c;
	vector<CubeVertex> cv;
	vector<Point3D> newtriangs;
	//#pragma omp parallel for
	for (int z = 0; z < depth; z++) {
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				cv.resize(0); //reset the vertices

				cv.push_back(getVertexAt(x, y+1, z+1)); //Corner 0, according to Bourke's picture
				cv.push_back(getVertexAt(x+1, y+1, z+1));
				cv.push_back(getVertexAt(x+1, y+1, z));
				cv.push_back(getVertexAt(x, y+1, z));
				cv.push_back(getVertexAt(x, y, z+1));
				cv.push_back(getVertexAt(x+1, y, z+1));
				cv.push_back(getVertexAt(x+1, y, z));
				cv.push_back(getVertexAt(x, y, z));
				c.setVertices(cv);
				c.getTriangles(0.5, triangles, normals, threshold);
			}
		}
	}
}
