#include "SpatialGrid.h"


SpatialGrid::SpatialGrid(int s, float h)
{
	this->sideLength = s * h;
	this->boxLength = h;
	this->start = Point3D(-sideLength/2, -sideLength/2, -sideLength/2); //the start of the boxes (corner of grid[0][0][0])
	this->numEdgeBoxes = s;
	//is at the left bottom far corner of the grid

	//x axis size
	grid.resize(s);
	for(int i=0;i<s;i++) 
	{
		//y axis size
		grid[i].resize(s);
			for(int j=0;j<s;j++)
			{
				//z axis size
				grid[i][j].resize(s);
			}
	}
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
	vector<Particle> list;

	//oh god
	while(i <= 1){
		int j=-1;
		while(j <= 1){
			int k=-1;
			while(k <= 1){
				if(!(xindex+i < 0 || xindex+i > this->numEdgeBoxes || yindex+j < 0 || yindex+j > this->numEdgeBoxes 
					|| zindex + k < 0 || zindex + k > this->numEdgeBoxes))
				{
					vector<Particle> thisBox = grid[xindex+i][yindex+j][zindex+k];
					int l = 0;
					while(l < thisBox.size()){
						list.push_back(thisBox[l]);
						l++;
					}
				}
				k++;
			}
			j++;
		}
		i++;
	}

	return list;

}

void SpatialGrid::updateBoxes(){
	int i = 0;
	while (i < this->numEdgeBoxes){
		int j = 0;
		while(j < this->numEdgeBoxes){
			int k = 0;
			while(k < this->numEdgeBoxes){
				int l = 0;
				while(l < grid[i][j][k].size()){
					Particle p = grid[i][j][k][l];
					Point3D pos = p.getPosition();
					float x = pos.getX();
					float y = pos.getY();
					float z = pos.getZ();
					if(abs(x) > sideLength || abs(y) > sideLength || abs(z) > sideLength){//moved particle is outside of spatial grid
						//cant add the particle!
						//uhhh not sure what to do here
					}
					int xindex = (int)floor((x+this->start.getX())/this->numEdgeBoxes);//need to check if this is right
					int yindex = (int)floor((y+this->start.getY())/this->numEdgeBoxes);
					int zindex = (int)floor((z+this->start.getZ())/this->numEdgeBoxes);
					if(xindex != i || yindex != j || zindex != k){
						grid[i][j][k].erase(grid[i][j][k].begin() + l);
						grid[xindex][yindex][zindex].push_back(p);
					}
					l++;
				}
				k++;
			}
			j++;
		}
		i++;
	}
}