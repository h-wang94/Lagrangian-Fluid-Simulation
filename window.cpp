/*
 * window.cpp
 *
 *  Created on: Dec 2, 2013
 *      Author: Owner
 */
#include "Particle.h"
#include "Vector.h"
#include "Point3D.h"
#include "ParticleSystem.h"
#include "SpatialGrid.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>

#include <GL/glut.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif
using namespace std;

std::vector<Point3D> vertexes, vt;
std::vector<Vector>normals;
std::vector< std::vector<int> > vertexIndexes;
ParticleSystem pSystem;
float currentTime = 0;
float dt;
float totalTime;

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

}

void testParticleSystem(unsigned const int numParticles, const float dt, const float totalTime, const bool objFile) {
  Particle p1 = Particle();
  unsigned int i;
  srand(time(NULL)); // for rand function
  if (objFile) {
    for(i = 0; i < numParticles; i++) {
      p1 = Particle(0.02, 0.0, 3.0, 998.29, 0, 3.5, vertexes[i], Vector(0, 0, 0));
      pSystem.addParticle(p1);
    }
  } else {
    for(i = 0; i < numParticles; i++) {
      p1 = Particle();
      pSystem.addParticle(p1);
    }
  }
  pSystem.initialize(dt);
  // dont need this now that we have opengl.
  // but we can use this if we don't want to display...i guess
  /*for(float j = 0; j < totalTime; j+=dt) {*/
    //cout << "//===========================================//" << endl
      //<< "// totalTime: " << j << "                                 //" << endl
      //<< "//===========================================//" << endl;
    //pSystem.update(dt);
  /*}*/
}

//==============================================================================
// Count the number of slashes to determine format of obj file inputs for faces
// =============================================================================
unsigned int countSlashes(std::string& s) {
  unsigned int counter = 0;
  unsigned int position = 0;
  for(position = 0; position < s.length(); position++) {
    if(s[position] == '/') {
      counter++;
    }
  }
  return counter;
}

//==============================================================================
// Parses the obj file.
// =============================================================================
void parseObjLine(std::vector<std::string>& splitline, std::vector<int>& temp) {
  float x, y, z;
  if(splitline.size() == 0) {
  }
  // reading vertex
  else if(!splitline[0].compare("v")) {
    if(splitline.size() == 4) {
      x = atof(splitline[1].c_str());
      y = atof(splitline[2].c_str());
      z = atof(splitline[3].c_str());
      Point3D vertex = Point3D(x, y, z);
      vertexes.push_back(vertex);
    }
    else {
      cout << "Error in number of parameters for V. Number of inputs should be 3. Has " << splitline.size() - 1 << endl;
    }
  }
  // don't have anything that will really care about this yet
  else if(!splitline[0].compare("vt")) {
    x = atof(splitline[1].c_str());
    y = atof(splitline[2].c_str());
    if(splitline.size() == 3) {
      z = 0;
    } else {
      z = atof(splitline[3].c_str());
    }
    Point3D vertex = Point3D(x, y, z);
    vt.push_back(vertex);
  }
  // vertex normals
  else if(!splitline[0].compare("vn")) {
    if(splitline.size() == 4) {
      x = atof(splitline[1].c_str());
      y = atof(splitline[2].c_str());
      z = atof(splitline[3].c_str());
      Vector normal = Vector(x, y, z);
      normals.push_back(normal);
    }
    else {
      cout << "Error in number of parameters for V. Number of inputs should be 3. Has " << splitline.size() - 1 << endl;
    }
  }
  // defines which vertexes will form a face
  else if(!splitline[0].compare("f")) { // will get the indexes of the vertexes, vt, and normals needed
    unsigned int numPoints = splitline.size() - 1;
    unsigned int numSlashes = countSlashes(splitline[1]);
    // correspond to the indexes
    // ==============================
    // Not using normals
    // ==============================
    //unsigned normalA, normalB, normalC, normal
    unsigned pointA, pointB, pointC, pointD;
    unsigned int position;
    if (numSlashes == 0) { // in the form "f v v v" with possibility of 4th
      temp.clear();
      pointA = atoi(splitline[1].c_str());
      pointB = atoi(splitline[2].c_str());
      pointC = atoi(splitline[3].c_str());
      temp.push_back(pointA);
      temp.push_back(pointB);
      temp.push_back(pointC);
      if (numPoints == 4) {
        pointD = atoi(splitline[4].c_str());
        temp.push_back(pointD);
      }
      vertexIndexes.push_back(temp);

    }
    else if (numSlashes == 1) { // in the form "f v/vt v/vt v/vt" with possibilty of 4th
      temp.clear();
      position = splitline[1].find("/");
      pointA = atoi(splitline[1].substr(0, position).c_str());
      position = splitline[2].find("/");
      pointB = atoi(splitline[2].substr(0, position).c_str());
      position = splitline[3].find("/");
      pointC = atoi(splitline[3].substr(0, position).c_str());
      temp.push_back(pointA);
      temp.push_back(pointB);
      temp.push_back(pointC);
      if (numPoints == 4) { // sometimes has 4 points
        position = splitline[4].find("/");
        pointD = atoi(splitline[4].c_str());
        temp.push_back(pointD);
      }
      vertexIndexes.push_back(temp);
    }

    else if (numSlashes == 2) { // in the form "f v/vt/n v/vt/n v/vt/n" or "f v//n v//n v//n" with possibility of 4th
      if(splitline[1].find("//")) {  // need to store normal
        temp.clear();
        position = splitline[1].find("/");
        pointA = atoi(splitline[1].substr(0, position).c_str());
        //normalA = atoi(splitline[1].substr(position + 2).c_str());
        position = splitline[2].find("/");
        pointB = atoi(splitline[2].substr(0, position).c_str());
        //normalB = atoi(splitline[2].substr(position + 2).c_str());
        position = splitline[3].find("/");
        pointC = atoi(splitline[3].substr(0, position).c_str());
        //normalC = atoi(splitline[3].substr(position + 2).c_str());
        temp.push_back(pointA);
        temp.push_back(pointB);
        temp.push_back(pointC);
        if (numPoints == 4) { // sometimes has 4 points
          position = splitline[4].find("/");
          pointD = atoi(splitline[4].c_str());
          //normalD = atoi(splitline[4].substr(position + 2).c_str());
          temp.push_back(pointD);
        }
        vertexIndexes.push_back(temp);
      }
      else { // in the form "f v//n v//n v//n" with possibility of 4th
        temp.clear();
        position = splitline[1].find("/");
        pointA = atoi(splitline[1].substr(0, position).c_str());
        position = splitline[1].find("/", position + 1);
        //normalA = atoi(splitline[1].substr(position + 1).c_str())

        position = splitline[2].find("/");
        pointB = atoi(splitline[2].substr(0, position).c_str());
        position = splitline[2].find("/", position + 1);
        //normalB = atoi(splitline[2].substr(position + 1).c_str());

        position = splitline[3].find("/");
        pointC = atoi(splitline[3].substr(0, position).c_str());
        position = splitline[3].find("/", position + 1);
        //normalC = atoi(splitline[3].substr(position + 1).c_str());
        temp.push_back(pointA);
        temp.push_back(pointB);
        temp.push_back(pointC);
        if (numPoints == 4) { // sometimes has 4 points
          position = splitline[4].find("/");
          pointD = atoi(splitline[4].substr(0, position).c_str());
          position = splitline[4].find("/", position + 1);
          //normalD = atoi(splitline[4].substr(position + 1).c_str());
          temp.push_back(pointD);
        }
        vertexIndexes.push_back(temp);
      }
    }
  }
}


void readInput(std::string fileName) {
  if (!(fileName.find(".obj") != std::string::npos)) {
    cout << "Error! Not an .OBJ file" << endl;
    exit(1);
  }
  // inputted file is a .obj file
  // start reading file
  std::ifstream inpfile(fileName.c_str());
  if (!inpfile.is_open()) {
    cout << "Error in opening file" << endl;
    exit(1);
  } else {
    cout << "//===================================//" << endl;
    cout << "// Begin reading " << fileName.c_str() << endl;
    std::string line;
    std::vector<int> temp;
    while(inpfile.good()) {
      std::vector<std::string> splitline;
      std::string buf;
      std::getline(inpfile,line);
      std::stringstream ss(line);
      while (ss >> buf) {
        splitline.push_back(buf);
      }
      // for obj files
      parseObjLine(splitline, temp);
    }
  }
  cout << "// Done reading " << fileName.c_str() << endl
    << "//===================================//" << endl
    << "//===================================//" << endl
    << "// Basic Info                        //" << endl
    << "//===================================//" << endl;
  cout << "// File name: " << fileName << endl;
  cout << "// Number of vertices: " << vertexes.size() << endl;
  cout << "//===================================//" << endl;
}

void determineFunction(int argc, char *argv[]) {
  bool objFile = false;
  unsigned int numParticles = 0;
  dt = 0;
  totalTime = 0;;
  std::string fileName;
  if (argc == 1) {
    numParticles = 1;
    dt = 0.001;
    totalTime = 1;
  }
  else if (argc == 4) {
    numParticles = std::atoi(argv[1]);
    dt = std::atof(argv[2]);
    totalTime = std::atof(argv[3]);
  }
  else if (argc == 2) {
    fileName = std::string(argv[1]);
    readInput(fileName); 
    objFile = true;
    numParticles = vertexes.size();
    dt = 0.001;
    totalTime = 1;
  }
  else if (argc == 5) {
    fileName = std::string(argv[1]);
    readInput(fileName);
    objFile = true;
    numParticles = vertexes.size();
    dt = std::atof(argv[2]);
    totalTime = std::atof(argv[3]); 
  }
  else {
    cout << "Incorrect number of parameters!" << endl;
    exit(1);
  }
  testParticleSystem(numParticles, dt, totalTime, objFile);
}


/* Display is updated/rendered here. */
void displayFunc() {
	glClear(GL_COLOR_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// rendering stuff goes here
	//for each particle
	glBegin(GL_POINTS);
	glColor4f(0.8, 0.2, 0.0, 0.5);
  Particle* p;
	for (unsigned int i = 0; i < pSystem.getParticles().size(); i+=1) {
    p = pSystem.getParticle(i);
    glVertex3f(p->getPosition().getX(), p->getPosition().getY(), p->getPosition().getZ());
	}
  if (currentTime < totalTime) {
    pSystem.update(dt);
    currentTime += dt;
  }
	glEnd();
	glFlush();
	glutSwapBuffers();
}

void keyboard(unsigned char key, int x, int y) {
  switch(key) {
    case ' ': // allow spacebar to end the program
      exit(0);
      break;
    default:
      break;
  }
  glutPostRedisplay();
}
int main(int argc, char** argv) {
  pSystem = ParticleSystem(Vector(0,-9.8, 0));
  determineFunction(argc, argv);

  glutInit(&argc, argv);
  glutInitWindowSize(600, 600);
  glutInitWindowPosition(0, 0);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
  glutCreateWindow("CS184 Final");
  glClearColor(0.0, 0.0, 0.0, 1.0);
  glEnable(GL_POINT_SMOOTH);
  glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);
  glPointSize(5.0f);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_BLEND);
  glutKeyboardFunc(keyboard);

  glutDisplayFunc(displayFunc);
  glutIdleFunc(displayFunc); //if we do user interaction
  glutMainLoop();
  return 0;
}


