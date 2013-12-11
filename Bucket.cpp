#include "Bucket.h"


Bucket::Bucket(void)
{
}

Bucket::Bucket(Point3D pos)
{
	this->position = pos;
}

Bucket::~Bucket(void)
{
}

std::vector<Particle> Bucket::getContents(){
	return this->contents;
}

std::vector<int> Bucket::getNeighbors(){
	return this->neighbors;
}

std::vector<Bucket*> Bucket::getBuckets(){
	return this->buckets;
}

Point3D Bucket::getPosition(){
	return this->position;
}

void Bucket::addParticle(Particle p){
	this->contents.push_back(p);
}

void Bucket::setContents(std::vector<Particle> list){
	this->contents = list;
}

void Bucket::setNeighbors(std::vector<int> list){
	this->neighbors = list;
}

void Bucket::addNeighbor(int n){
	this->neighbors.push_back(n);
}

void Bucket::removeParticle(Point3D position){
	int c = 0;
	while(c < this->contents.size()){
		//if(this->contents[c].getPosition
		c++;
	}
}
