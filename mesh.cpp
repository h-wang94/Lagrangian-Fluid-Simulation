/*
 * mesh.cpp
 *
 *  Created on: Dec 11, 2013
 *      Author: Owner
 */

#include "mesh.h"

Mesh::Mesh() {
	width = 0;
	height = 0;
	depth = 0;
}

Mesh::Mesh(int width, int height, int depth, int step) {
	this->width = width;
	this->height = height;
	this->depth = depth;
	mesh.resize(width * height * depth); //ow
	for (int z = 0; z < depth; z++) {
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				setVertexAt(x, y, z, CubeVertex(Water(Point3D(x * step, y * step, z * step))));
			}
		}
	}
}

CubeVertex Mesh::getVertexAt(int x, int y, int z) {
	return mesh[(z * width * height) + (y * width) + x]; //wat
}

void Mesh::setVertexAt(int x, int y, int z, const CubeVertex &c) {
	mesh[(z * width * height) + (y * height) + x] = c;
}

void Mesh::updateColors() {
	for (int z = 0; z < depth; z++) {
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				CubeVertex cv = getVertexAt(x, y, z);
				Particle p = cv.getParticle();
				cv.setColor(pSystem.colorFunction(p));
			}
		}
	}
}

vector<vector<Point3D> > Mesh::marchingCubes() {
	Cube c;
	vector<CubeVertex> cv;
	vector<vector<Point3D> > triangles;
	for (int z = 0; z < depth - 1; z++) {
		for (int y = 0; y < height - 1; y++) {
			for (int x = 0; x < width - 1; x++) {
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
				triangles.push_back(c.getTriangles(0.5));
			}
		}
	}
	return triangles;
}
