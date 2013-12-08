#ifndef _SpatialGrid_h_
#define _SpatialGrid_h_
#include <vector>
#include "Particle.h"

class SpatialGrid
{
public:
	SpatialGrid(const int s, const float h);//establish grid with h box lengths of s*h side lengths of a huge cube grid thingy at point 0,0,0
	SpatialGrid();
	void addParticle(Particle p); //add particle to list of particles in box
	void updateBoxes(std::vector<Particle> particles); //when particles move, they might change boxes, so we need to update boxes
	std::vector<Particle> getNeighbors(Particle p); //get the neighbor particles in the box and in the surrounding boxes
	
private:
	std::vector<std::vector<std::vector<std::vector<Particle> > > > grid; //the grid consists of 3d boxes with a list of particles in each box
	float sideLength;
	float boxLength;
	int numEdgeBoxes;
	Point3D start;

};

#endif

