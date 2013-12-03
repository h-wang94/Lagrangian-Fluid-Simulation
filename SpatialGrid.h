#ifndef _SpatialGrid_h_
#define _SpatialGrid_h_
#include <vector>
#include "Particle.h"

class SpatialGrid
{
public:
	SpatialGrid(int s, float h);//establish grid with h box lengths of s*x side lengths of a huge cube grid thingy at point 0,0,0
	void addParticle(Particle p); //add particle to list of particles in box
	void updateBoxes(); //when particles move, they might change boxes, so we need to update boxes
	std::vector<Particle> getNeighbors(Particle p); //get the neighbor particles in the box and in the surrounding boxes
	
private:
	std::vector<std::vector<std::vector<std::vector<Particle> > > > grid;
	float sideLength;
	float boxLength;
	int numEdgeBoxes;
	Point3D start;

};

#endif
