#ifndef _SpatialGrid_h_
#define _SpatialGrid_h_
#include <vector>
#include "Particle.h"
#include <map>
#include "Bucket.h"

class SpatialGrid
{
public:
	SpatialGrid(const float h);//establish grid with h box lengths
	SpatialGrid();
	int addParticle(Particle p); //add particle to list of particles in box
	std::vector<Particle> updateBoxes(std::vector<Particle> particles); //when particles move, they might change boxes, so we need to update boxes
	std::vector<Particle> getNeighbors(Particle p); //get the neighbor particles in the box and in the surrounding boxes
	void clear();
	
private:
	//std::vector<std::vector<std::vector<std::vector<Particle> > > > grid; //the grid consists of 3d boxes with a list of particles in each box
	std::map<int, Bucket> hashMap;
	float sideLength;
	float boxLength;
	int numEdgeBoxes;
	int p1;
	int p2;
	int p3;
	Point3D start;
	bool newPositionsSet;

};

#endif

