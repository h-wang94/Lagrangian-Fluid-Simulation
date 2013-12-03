#include "SpatialGrid.h"
#include <cmath>


SpatialGrid::SpatialGrid(int s, float h)
{
	this->sideLength = s * h;
	this->boxLength = h;
	this->start = Point3D(-sideLength/2, -sideLength/2, -sideLength/2); //the start of the boxes (corner of grid[0][0][0])
	this->numEdgeBoxes = s;
	//is at the left bottom far corner of the grid
}

void SpatialGrid::addParticle(Particle p){
	Point3D pos = p.getPosition();
	float x = pos.getX();
	float y = pos.getY();
	float z = pos.getZ();
	if(abs(x) > sideLength || abs(y) > sideLength || abs(z) > sideLength){//added particle is outside of spatial grid
		//cant add the particle!
		return;
	}
	int xindex = (int)floor((x+this->start.getX())/this->numEdgeBoxes);//need to check if this is right
	int yindex = (int)floor((y+this->start.getY())/this->numEdgeBoxes);
	int zindex = (int)floor((z+this->start.getZ())/this->numEdgeBoxes);
	grid[xindex][yindex][zindex].push_back(p);
}

std::vector<Particle> SpatialGrid::getNeighbors(Particle p){
	Point3D pos = p.getPosition();
	float x = pos.getX();
	float y = pos.getY();
	float z = pos.getZ();
	int xindex = (int)floor((x+this->start.getX())/this->numEdgeBoxes);//need to check if this is right
	int yindex = (int)floor((y+this->start.getY())/this->numEdgeBoxes);
	int zindex = (int)floor((z+this->start.getZ())/this->numEdgeBoxes);

	int i = -1;
	int j = -1;
	int k = -1;
	vector<Particle> list;

	//oh god
	while(i <= 1){
		while(j <= 1){
			while(k <= 1){
				vector<Particle> thisBox = grid[xindex+i][yindex+j][zindex+k];
				unsigned int l = 0;
				while(l < thisBox.size()){
					list.push_back(thisBox[l]);
					l++;
				}
				l = 0;
				k++;
			}
			j++;
		}
		i++;
	}

}

void SpatialGrid::updateBoxes(){

}
