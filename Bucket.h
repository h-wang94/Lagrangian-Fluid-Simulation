#ifndef _Bucket_h_
#define _Bucket_h_
#include <vector>
#include "Particle.h"
class Bucket // A bucket stores the list of particles in the grid box and the buckets neighbors
{
public:
	Bucket(void);
	Bucket(Point3D pos);
	~Bucket(void);
	std::vector<Particle> getContents();
	void setContents(std::vector<Particle> list);
	std::vector<int> getNeighbors();
	std::vector<Bucket*> getBuckets();
	Point3D getPosition();
	void addParticle(Particle p);
	void addNeighbor(int n);
	void setNeighbors(std::vector<int> list);
	void removeParticle(Point3D position);

private:
	std::vector<Particle> contents;
	std::vector<int> neighbors;
	std::vector<Bucket*> buckets;
	Point3D position;
};
#endif

