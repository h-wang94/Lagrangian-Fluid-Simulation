#include "SpatialGrid.h"
#include <cmath>


SpatialGrid::SpatialGrid(int s, float h)//s is the number of boxes for a side, h is the smoothing radius of our weighting function
{
	this->sideLength = s * h;
	this->boxLength = h;
	this->start = Point3D(-s * h/2, -s * h/2, -s * h/2); //the start of the boxes (corner of grid[0][0][0])
	this->numEdgeBoxes = s;
	//cout<<start<<endl;
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

SpatialGrid::SpatialGrid(){

}

void SpatialGrid::addParticle(Particle p){
	Point3D pos = p.getPosition();
	float x = pos.getX();
	float y = pos.getY();
	float z = pos.getZ();
	if(abs(x) > sideLength/2 || abs(y) > sideLength/2 || abs(z) > sideLength/2){//added particle is outside of spatial grid
		cout<<"Particle outside of boundary!"<<endl;
		//cant add the particle!
		return;
	}
	int xindex = (int)floor((x-this->start.getX())/(this->boxLength));
	int yindex = (int)floor((y-this->start.getY())/(this->boxLength));
	int zindex = (int)floor((z-this->start.getZ())/(this->boxLength));
	/*cout<<xindex<<", "<<yindex<<", "<<zindex<<"."<<endl;
	cout<<pos<<endl;
	cout<<sideLength<<endl;
	cout<<start<<endl;*/
	grid[xindex][yindex][zindex].push_back(p);
}

std::vector<Particle> SpatialGrid::getNeighbors(Particle p){
	Point3D pos = p.getPosition();
	float x = pos.getX();
	float y = pos.getY();
	float z = pos.getZ();
	int xindex = (int)floor((x-this->start.getX())/(this->boxLength));
	int yindex = (int)floor((y-this->start.getY())/(this->boxLength));
	int zindex = (int)floor((z-this->start.getZ())/(this->boxLength));

	int i = -1;
	vector<Particle> list;

	//oh god
	while(i <= 1){
		int j=-1;
		while(j <= 1){
			int k=-1;
			while(k <= 1){
				//cout<<xindex<<", "<<yindex<<", "<<zindex<<endl;
				if(!(xindex+i < 0 || xindex+i >= this->numEdgeBoxes || yindex+j < 0 || yindex+j >= this->numEdgeBoxes 
					|| zindex + k < 0 || zindex + k >= this->numEdgeBoxes))
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

void SpatialGrid::updateBoxes(std::vector<Particle> particles){
	int n = 0;
	while(n < particles.size()){
		Particle p = particles[n];
		//cout<<p<<endl;
		Point3D pos = p.getPosition();
		float x = pos.getX();
		float y = pos.getY();
		float z = pos.getZ();
		float oldX = p.getOldPosition().getX();
		float oldY = p.getOldPosition().getY();
		float oldZ = p.getOldPosition().getZ();

		int xindex = (int)floor((x-this->start.getX())/(this->boxLength));
		int yindex = (int)floor((y-this->start.getY())/(this->boxLength));
		int zindex = (int)floor((z-this->start.getZ())/(this->boxLength));
		int i = (int)floor((oldX-this->start.getX())/(this->boxLength));
		int j = (int)floor((oldY-this->start.getY())/(this->boxLength));
		int k = (int)floor((oldZ-this->start.getZ())/(this->boxLength));

		/*cout<<"--------"<<endl;
		cout<<i<<", "<<j<<", "<<k<<endl;
		cout<<xindex<<", "<<yindex<<", "<<zindex<<endl;
		cout<<"--------"<<endl;*/

		//cout<<p<<endl;
		//cout<<p.getOldPosition()<<endl;
		if(abs(x) > sideLength/2 || abs(y) > sideLength/2 || abs(z) > sideLength/2){//moved particle is outside of spatial grid
			cout<<"Particle outside of boundary!"<<endl;
			//particle is OOB, bye particle
			//grid[i][j][k].erase(grid[i][j][k].begin() + l);
		}
		else{

			if(xindex != i || yindex != j || zindex != k){
				//cout<<"omg"<<endl;
				int l = 0;
					while (l < grid[i][j][k].size()){
						if(grid[i][j][k][l].getPosition().getX() == oldX 
							&& grid[i][j][k][l].getPosition().getY() == oldY
							&& grid[i][j][k][l].getPosition().getZ() == oldZ){
							grid[i][j][k].erase(grid[i][j][k].begin() + l);
						}
						l++;
					grid[xindex][yindex][zindex].push_back(p);
				}
			}
		}
		n++;
	}

	/*while (i < this->numEdgeBoxes){
		int j = 0;
		while(j < this->numEdgeBoxes){
			int k = 0;
			while(k < this->numEdgeBoxes){
				int l = 0;
				while(l < oldGrid.grid[i][j][k].size()){
					Particle p = oldGrid.grid[i][j][k][l];
					//cout<<p<<endl;
					Point3D pos = p.getPosition();
					float x = pos.getX();
					float y = pos.getY();
					float z = pos.getZ();
					if(abs(x) > sideLength || abs(y) > sideLength || abs(z) > sideLength){//moved particle is outside of spatial grid
						//particle is OOB, bye particle
						//grid[i][j][k].erase(grid[i][j][k].begin() + l);
					}
					else{
						int xindex = (int)floor((x-this->start.getX())/(this->boxLength));
						int yindex = (int)floor((y-this->start.getY())/(this->boxLength));
						int zindex = (int)floor((z-this->start.getZ())/(this->boxLength));
						if(xindex != i || yindex != j || zindex != k){
							//cout<<"omg"<<endl;
							grid[i][j][k].erase(grid[i][j][k].begin() + l);
							grid[xindex][yindex][zindex].push_back(p);
						}
					}
					l++;
				}
				k++;
			}
			j++;
		}
		i++;
	}*/
}

