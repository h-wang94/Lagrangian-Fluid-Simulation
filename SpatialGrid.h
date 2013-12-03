#ifndef _SpatialGrid_h_
#define _SpatialGrid_h_
#include <vector>
#include "Particle.h"

class SpatialGrid
{
public:
	SpatialGrid(float s, float h);//establish grid with h box lengths and s length sides of a huge cube grid thingy at point 0,0,0
	void addParticle(Point3D p); //add particle to list of particles in box
	
private:
	std::vector<std::vector<std::vector<std::vector<Particle> > > > grid;

};

#endif

