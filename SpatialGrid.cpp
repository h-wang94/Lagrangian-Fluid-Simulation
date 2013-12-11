#include "SpatialGrid.h"
#include <omp.h>
#include <cmath>
#include <cstdlib>


SpatialGrid::SpatialGrid(const float h)//s is the number of boxes for a side, h is the smoothing radius of our weighting function
{
	this->boxLength = h;

	this->p1 = 10619863;
	this->p2 = 94418953;
	this->p3 = 54018521;

	this->newPositionsSet = false;
}

SpatialGrid::SpatialGrid(){

}

void SpatialGrid::clear(){

}

int SpatialGrid::addParticle(Particle p){
	Point3D pos = p.getPosition();
	float x = pos.getX();
	float y = pos.getY();
	float z = pos.getZ();

	int xindex = (int)floor((x+10)/(this->boxLength));
	int yindex = (int)floor((y+10)/(this->boxLength));
	int zindex = (int)floor((z+10)/(this->boxLength));

	int thisid = ((xindex*p1) ^ (yindex*p2) ^ (zindex*p3)); //no idea where this particle is in the grid yet

	if(this->hashMap.find(thisid) != this->hashMap.end()){//such a bucket exists
		std::map<int,Bucket>::iterator it;
		it = this->hashMap.find(thisid);
		it->second.addParticle(p);
	}
	else{//the bucket does not exist
		Bucket b =Bucket(Point3D(xindex, yindex, zindex));
		b.addParticle(p);
		hashMap.insert ( std::pair<int, Bucket>(thisid,b) );
	}

	std::map<int,Bucket>::iterator it;
	it = this->hashMap.find(thisid);
	std::vector<int> ids = it->second.getNeighbors();

	//either way we need to check neighbors of the bucket if they don't exist
	if(it->second.getNeighbors().size() != 26){//26 surrounding boxes total
		int i = -1;
		while(i <= 1){
			int j = -1;
			while(j <= 1){
				int k = -1;
				while(k <= 1){
					int xindex2 = xindex + i;
					int yindex2 = yindex + j;
					int zindex2 = zindex + k;
					int id = ((xindex2*p1) ^ (yindex2*p2) ^ (zindex2*p3)); //instantiate neighbor ids
					//if(this->hashMap.find(id) != this->hashMap.end()){//such a bucket exists
					//}
					//else{//the bucket does not exist
						if (!(i==0 && j==0 && k==0)){
							//this->hashMap[id] = Bucket(Point3D(xindex2, yindex2, zindex2));
							ids.push_back(id);
						}
					//}
					k++;
				}
				j++;
			}
			i++;
		}
	}

	return thisid;

}

std::vector<Particle> SpatialGrid::getNeighbors(Particle p){
	int thisid = p.getHashID();//((xindex*p1) ^ (yindex*p2) ^ (zindex*p3));
	vector<Particle> list;
	Bucket it = this->hashMap.find(thisid)->second;

	//if(it != this->hashMap.end()){//such a bucket exists

		std::vector<int> neighbors = it.getNeighbors();
		std::vector<Particle> parts = it.getContents();
		unsigned int i = 0;
//#pragma omp parallel for firstprivate(i)
		for(i = 0; i < parts.size(); i++){
			/*if(!(it->second.getContents()[i].getPosition().getX() == x &&
							it->second.getContents()[i].getPosition().getY() == y &&
							it->second.getContents()[i].getPosition().getZ() == z)){*/
				list.push_back(parts[i]);
			//}
		}

		unsigned int j = 0;
//#pragma omp parallel for firstprivate(j)
		for(j = 0; j < neighbors.size(); j++){
			unsigned int k = 0;
			//std::map<int,Bucket>::iterator it2 = this->hashMap.find(neighbors[j]);
      if(this->hashMap.find(neighbors[j]) != this->hashMap.end()) {
        parts = this->hashMap.find(neighbors[j])->second.getContents();
        for(k=0; k < parts.size(); k++){
          list.push_back(parts[k]);
        }
      }
		}
	//}
	//else{//the bucket does not exist
		// it should never go in here since particles are in the grid/bucket in the first place
		/*cout<<"ERROR: Checking neighbors of particle not in grid"<<endl;
		int breakpoint = 0;
		cin >> breakpoint;
		exit(666);*/
	//}
	//cout << "DONE WITH GET NEIGHBORS" << endl;
	return list;

}

std::vector<Particle> SpatialGrid::updateBoxes(std::vector<Particle> particles){
	
	if(newPositionsSet == false){//haven't even initialized the new positions yet
		unsigned int c = 0;
		while(c < particles.size()){
			Point3D pos = particles[c].getPosition();//particle now has a new position
			float x = pos.getX();
			float y = pos.getY();
			float z = pos.getZ();
			Point3D OLDpos = particles[c].getOldPosition();
			float OLDx = OLDpos.getX();
			float OLDy = OLDpos.getY();
			float OLDz = OLDpos.getZ();
			int xindex = (int)floor((x+10)/(this->boxLength));//these are now new position indices
			int yindex = (int)floor((y+10)/(this->boxLength));
			int zindex = (int)floor((z+10)/(this->boxLength));
			int thisid = particles[c].getHashID();//((OLDxindex*p1) ^ (OLDyindex*p2) ^ (OLDzindex*p3));//this is the OLD HASH ID
			int newid = ((xindex*p1) ^ (yindex*p2) ^ (zindex*p3));//this is the NEW HASH ID we need to set
			particles[c].setHashID(newid);
			
			/*cout<<"the list parts hash id: "<<thisid<<endl;
			cout<<"the calculated new position hash id: "<<newid<<endl;
			cout<<"the calculated old position hash id: "<<((OLDxindex*p1) ^ (OLDyindex*p2) ^ (OLDzindex*p3))<<endl;
			cout<<"the old position: "<<endl<<OLDpos<<endl;
			cout<<"the new position: "<<endl<<pos<<endl;
			int breakpoint;
			cin >> breakpoint;*/

			if(thisid != newid){
				particles[c].setHashID(newid);
				//need to put the particle into the new bucket and erase the other particle
				//add particle in new position
				if(this->hashMap.find(newid) != this->hashMap.end()){//such a bucket exists
					std::map<int,Bucket>::iterator it;
					it = this->hashMap.find(newid);
					it->second.addParticle(particles[c]);
				}
				else{//the bucket does not exist
					Bucket b =Bucket(Point3D(xindex, yindex, zindex));
					b.addParticle(particles[c]);
					hashMap.insert ( std::pair<int, Bucket>(newid,b) );
				}

				std::map<int,Bucket>::iterator it;
				it = this->hashMap.find(newid);
				std::vector<int> ids = it->second.getNeighbors();
				//either way we need to check neighbors of the bucket if they don't exist
				if(it->second.getNeighbors().size() != 26){//26 surrounding boxes total
					int i = -1;
					while(i <= 1){
						int j = -1;
						while(j <= 1){
							int k = -1;
							while(k <= 1){
								int xindex2 = xindex + i;
								int yindex2 = yindex + j;
								int zindex2 = zindex + k;
								int id = ((xindex2*p1) ^ (yindex2*p2) ^ (zindex2*p3)); //instantiate neighbor ids
								//if(this->hashMap.find(id) != this->hashMap.end()){//such a bucket exists
									//do nothing
								//}
								//else{//the bucket does not exist
									if (!(i==0 && j==0 && k==0)){
										//this->hashMap[id] = Bucket(Point3D(xindex2, yindex2, zindex2));
										ids.push_back(id);
									}
								//}
								k++;
							}
							j++;
						}
						i++;
					}
					it->second.setNeighbors(ids);
				}


				//delete old particle
				if(this->hashMap.find(thisid) != this->hashMap.end()){//such a bucket exists
					std::map<int,Bucket>::iterator it;
					it = this->hashMap.find(thisid);
					std::vector<Particle> list = it->second.getContents();
					unsigned int d = 0;
					while(d < list.size()){
						if(list[d].getPosition().getX() == OLDx &&
							list[d].getPosition().getY() == OLDy &&
							list[d].getPosition().getZ() == OLDz){
								list.erase(list.begin()+d);
						}
						d++;
					}
					/*it->second.setContents(list);
					//update the list oldposition
					it = this->hashMap.find(newid);
					list = it->second.getContents();*/
					d = 0;
					while(d < list.size()){
						if(list[d].getOldPosition().getX() == OLDx &&
							list[d].getOldPosition().getY() == OLDy &&
							list[d].getOldPosition().getZ() == OLDz){
								list[d].setOldPosition(Point3D(x,y,z));
							}
					
						d++;
					}
					it->second.setContents(list);

				}
				else{//the bucket does not exist
					//uhHHH?????? there should be an old particle so it shouldn't go in here
					cout<<"ERROR: WEIRD CS ANOMALY"<<endl;
					int breakpoint = 0;
					cin >> breakpoint;
					exit(1000);
				}
			}
			else{
				//cout<<"NO BUCKET CHANGE"<<endl;
				//update the list
				//std::map<int,Bucket>::iterator it;
				//it = this->hashMap.find(thisid);
				std::vector<Particle> list = this->hashMap.find(thisid)->second.getContents();
				unsigned int d = 0;
				while(d < list.size()){
					if(list[d].getPosition().getX() == OLDx &&
						list[d].getPosition().getY() == OLDy &&
						list[d].getPosition().getZ() == OLDz){
							list[d]=particles[c];
						}
					
					d++;
				}
				this->hashMap.find(thisid)->second.setContents(list);
			}

			c++;
			
		}
	}
	return particles;
}

