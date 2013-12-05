#include "Particle.h"
#include "Vector.h"
#include "Point3D.h"
#include "ParticleSystem.h"
#include "SpatialGrid.h"
#include <iostream>
#include <cstdlib>
using namespace std;

void testSpatialGrid() {
  /*Particle p1 = Particle(1, 1, 3.5, 1, Vector(0,0,0), Point3D(0,0,0));*/
  //Particle p2 = Particle(1, 1, 3.5, 1, Vector(0,0,0), Point3D(1,0,0));
  //Particle p3 = Particle(1, 1, 3.5, 1, Vector(0,0,0), Point3D(0,0,1));
  //Particle p4 = Particle(1, 1, 3.5, 1, Vector(0,0,0), Point3D(0,-9, -9));
  /*Particle p5 = Particle(1, 1, 3.5, 1, Vector(0,0,0), Point3D(11,11, 11));*/
  Particle p1 = Particle();
  Particle p2 = Particle();
  Particle p3 = Particle();
  Particle p4 = Particle();
  Particle p5 = Particle();
  SpatialGrid grid = SpatialGrid(10, 10);
  grid.addParticle(&p1);
  grid.addParticle(&p2);
  grid.addParticle(&p3);
  grid.addParticle(&p4);
  grid.addParticle(&p5);
  std::vector<Particle*> list = grid.getNeighbors(p1);
  unsigned int j = 0;
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
  cout<<p4<<endl;
  std::vector<Particle*> list2 = grid.getNeighbors(p1);
  j = 0;
  while (j < list2.size()){
    cout << *list2[j] << endl;
    j++;
  }

  /*unsigned int i;*/
  /*cin>>i;*/
}

void testParticleSystem(const int argc, char* argv[]) {
  unsigned int numParticles;
  float dt;
  float time;
  if(argc == 1) {
    numParticles = 1;
    dt = 0.1;
    time = 0.5;
  }
  else if (argc == 4) {
    numParticles = std::atoi(argv[1]);
    dt = std::atof(argv[2]);
    time = std::atof(argv[3]);
  }
  else {
    cout << "Error! Incorrect number of parameters" << endl;
    exit(1);
  }
  ParticleSystem system = ParticleSystem(Vector(0,0,-9.8));
  Particle p1 = Particle();
  //Particle p2 = Particle(2, 2, 1, 1, Vector(0,1,1), Point3D(1, 1, 1));
  for(unsigned int i = 0; i < numParticles; i++) {
    p1 = Particle();
    system.addParticle(p1);
  }
  system.initialize(dt);
  for(float j = 0; j < time; j+=dt) {
    cout << "//===========================================//" << endl
         << "// Time: " << j << "                                 //" << endl
         << "//===========================================//" << endl;
    system.update(dt);
  }
}

int main(int argc, char *argv[]) {
  cout << "//===========================================//" << endl
       << "// Fluid Simulation Project " << endl
       << "//===========================================//" << endl;
  
  testParticleSystem(argc, argv);
  //testSpatialGrid();
}


