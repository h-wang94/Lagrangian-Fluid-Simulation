#include "Particle.h"
#include "Vector.h"
#include "Point3D.h"
#include "ParticleSystem.h"
#include "SpatialGrid.h"
#include <iostream>
using namespace std;
int main() {
  cout << "//===========================================//" << endl
       << "// Fluid Simulation Project " << endl
       << "//===========================================//" << endl;

  // Testing print/debug statement
  /*Particle particle = Particle(1, 2, 3, Vector(4, 5, 6), Vector(7, 8, 9));*/
  /*cout << particle;*/

  /*Particle p1 = Particle(1, 1, 1, 1, Vector(0,0,0), Point3D(0,0,0));
  Particle p2 = Particle(1, 1, 1, 1, Vector(0,0,0), Point3D(1,0,0));
  Particle p3 = Particle(1, 1, 1, 1, Vector(0,0,0), Point3D(0,0,1));
  Particle p4 = Particle(1, 1, 1, 1, Vector(0,0,0), Point3D(0,-9, -9));
  Particle p5 = Particle(1, 1, 1, 1, Vector(0,0,0), Point3D(11,11, 11));
  SpatialGrid grid = SpatialGrid(10, 10);
  grid.addParticle(&p1);
  grid.addParticle(&p2);
  grid.addParticle(&p3);
  grid.addParticle(&p4);
  grid.addParticle(&p5);
  std::vector<Particle*> list = grid.getNeighbors(p1);
  int j = 0;
  while (j < list.size()){
	  cout << *list[j] << endl;
	  j++;
  }
  p4.setPosition(Point3D(10,10,10));
  cout<<"++++++++++++++++++++++++++++++++++++"<<endl;
  cout<<"++++++++++++++++++++++++++++++++++++"<<endl;
  cout<<"++++++++++++++++++++++++++++++++++++"<<endl<<endl;
  grid.updateBoxes();
  cout<<"++++++++++++++++++++++++++++++++++++"<<endl;
  cout<<"++++++++++++++++++++++++++++++++++++"<<endl;
  cout<<"++++++++++++++++++++++++++++++++++++"<<endl<<endl;
  //cout<<p4<<endl;
  std::vector<Particle*> list2 = grid.getNeighbors(p1);
  j = 0;
  while (j < list2.size()){
	  cout << *list2[j] << endl;
	  j++;
  }


  int i;
  cin>>i;*/



}
