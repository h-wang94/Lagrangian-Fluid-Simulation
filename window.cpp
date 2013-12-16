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
#include "GLToMovie.h"
#include "Film.h"
#include "mesh.h"
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


int frameWidth = 600, frameHeight = 600;
CGLToMovie recorder;
std::vector<Point3D> vertexes, vt;
std::vector<Vector>normals;
std::vector< std::vector<int> > vertexIndexes;
ParticleSystem pSystem;
Mesh mesh;
bool objFile = false;
bool record = false;
unsigned int numParticles;
float currentTime = 0;
float dt = 0;
float totalTime = 0;

bool bBoxScene = true;
bool bxScene = false;
bool bzScene = false;
bool bNoBound = false;

float maxX = 0.1f;
float maxY = 0.0f;
float maxZ = 0.1f;

bool useCubes = false;

void testSpatialGrid() {
  /*Particle p1 = Particle(1, 1, 3.5, 1, Vector(0,0,0), Point3D(0,0,0));*/
  //Particle p2 = Particle(1, 1, 3.5, 1, Vector(0,0,0), Point3D(1,0,0));
  //Particle p3 = Particle(1, 1, 3.5, 1, Vector(0,0,0), Point3D(0,0,1));
  //Particle p4 = Particle(1, 1, 3.5, 1, Vector(0,0,0), Point3D(0,-9, -9));
  /*Particle p5 = Particle(1, 1, 3.5, 1, Vector(0,0,0), Point3D(11,11, 11));*/
  std::vector<Particle> listofparts;
  Particle p1 = Particle();
  Particle p2 = Particle();
  Particle p3 = Particle();
  Particle p4 = Particle();
  Particle p5 = Particle();
  SpatialGrid grid = SpatialGrid(10);
  grid.addParticle(p1);
  grid.addParticle(p2);
  grid.addParticle(p3);
  grid.addParticle(p4);
  grid.addParticle(p5);
  listofparts.push_back(p1);
  listofparts.push_back(p2);
  listofparts.push_back(p3);
  listofparts.push_back(p4);
  listofparts.push_back(p5);
  std::vector<Particle> list = grid.getNeighbors(p1);
  unsigned int j = 0;
  while (j < list.size()){
    cout << list[j] << endl;
    j++;
  }
  p4.setPosition(Point3D(10,10,10));
  cout<<"++++++++++++++++++++++++++++++++++++"<<endl;
  cout<<"++++++++++++++++++++++++++++++++++++"<<endl;
  cout<<"++++++++++++++++++++++++++++++++++++"<<endl<<endl;
  int k = 0;
  while(k<listofparts.size()){
    listofparts[k] = grid.updateBoxes(listofparts[k]);
    k++;
  }
  cout<<"++++++++++++++++++++++++++++++++++++"<<endl;
  cout<<"++++++++++++++++++++++++++++++++++++"<<endl;
  cout<<"++++++++++++++++++++++++++++++++++++"<<endl<<endl;
  cout<<p4<<endl;
  std::vector<Particle> list2 = grid.getNeighbors(p1);
  j = 0;
  while (j < list2.size()){
    cout << list2[j] << endl;
    j++;
  }
  int breakpoint = 0;
  cin>>breakpoint;
}

void createParticles() {
  Particle p1 = Particle();
  unsigned int i;
  srand(time(NULL)); // for rand function
  if (objFile) {
    maxX = vertexes[0].getX();
    maxY = vertexes[0].getY();
    maxZ = vertexes[0].getZ();
    for(i = 0; i < vertexes.size(); i++) {
      //p1 = Particle(0.02, 0.0, 3.0, 998.29, 0, 3.5, vertexes[i], Vector(0, 0, 0));
      p1 = Water(vertexes[i]);
      if (maxX < vertexes[i].getX()) {
        maxX = vertexes[i].getX();
      }
      if (maxY < vertexes[i].getY()) {
        maxY = vertexes[i].getY();
      }
      if (maxZ < vertexes[i].getZ()) {
        maxZ = vertexes[i].getZ();
      }
      pSystem.addParticle(p1);
    }
    pSystem.setBoundaries(-maxX, -1, -maxZ, maxX + 1, maxY + 1, maxZ + 1);
  } else {
    for(i = 0; i < numParticles/2; i++) {
      p1 = Water();
      pSystem.addParticle(p1);
    }
    for(; i < numParticles; i++) {
      p1 = Mucus();
      pSystem.addParticle(p1);
    }
  }  
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


/*Parses commandline arguments. 
  Each optional flag should be followed by appropriate arguments for that flag. */
void determineFunction(int argc, char *argv[]) {
  numParticles = 1;
  dt = 0.001;
  totalTime = 1;
  for (int i = 1; i < argc;i++) {
    if (strcmp(argv[i], "-pNum") == 0) {
      i++;
      numParticles = std::atoi(argv[i]);
    } else if (strcmp(argv[i], "-dt") == 0) {
      i++;
      dt = std::atof(argv[i]);
    } else if (strcmp(argv[i], "-time") == 0) {
      i++;
      totalTime = std::atof(argv[i]);
    } else if (strcmp(argv[i], "-objInput") == 0) {
      i++;
      bNoBound = true;
      bBoxScene = false;
      objFile = true;
      readInput(string(argv[i]));
    } else if (strcmp(argv[i], "-width") == 0) {
      i++;
      frameWidth = std::atoi(argv[i]);
      recorder.setWidth(frameWidth);
    } else if (strcmp(argv[i], "-height") == 0) {
      i++;
      frameHeight = std::atoi(argv[i]);
      recorder.setHeight(frameHeight);
    } else if (strcmp(argv[i], "-record") == 0) {
      record = true;
    } else if (strcmp(argv[i], "-cubes") == 0) {
    	useCubes = true;
    } else {
      cerr << "Unknown argument: " << argv[i] << endl;
    }
  }
}

void boxScene() {
  glBegin(GL_LINE_STRIP);
  glColor3f(1.0, 1.0, 1.0);
  glVertex3f(-.11, -1, -.11);
  glVertex3f(-.11, 0, -.11);
  glVertex3f(-.11, 0, .11);
  glVertex3f(-.11, -1, .11);
  glVertex3f(-.11, -1, -.11);
  glVertex3f(.11, -1, -.11);
  glVertex3f(.11, 0, -.11);
  glVertex3f(-.11, 0, -.11);
  glVertex3f(-.11, -1, -.11);
  glEnd();
  glBegin(GL_LINE_STRIP);
  glColor3f(1.0, 1.0, 1.0);
  glVertex3f(.11, -1, .11);
  glVertex3f(-.11, -1, .11);
  glVertex3f(-.11, 0, .11);
  glVertex3f(.11, 0, .11);
  glVertex3f(.11, -1, .11);
  glVertex3f(.11, -1, -.11);
  glVertex3f(.11, 0, -.11);
  glVertex3f(.11, 0, .11);
  glVertex3f(.11, -1, .11);
  glEnd();
}

void xScene() {
  glBegin(GL_LINE_STRIP);
  glColor3f(1.0, 1.0, 1.0);
  glVertex3f(-1.01, -1, -.11);
  glVertex3f(-1.01, 0, -.11);
  glVertex3f(-1.01, 0, .11);
  glVertex3f(-1.01, -1, .11);
  glVertex3f(-1.01, -1, -.11);
  glVertex3f(.11, -1, -.11);
  glVertex3f(.11, 0, -.11);
  glVertex3f(-1.01, 0, -.11);
  glVertex3f(-1.01, -1, -.11);
  glEnd();
  glBegin(GL_LINE_STRIP);
  glColor3f(1.0, 1.0, 1.0);
  glVertex3f(.11, -1, .11);
  glVertex3f(-1.01, -1, .11);
  glVertex3f(-1.01, 0, .11);
  glVertex3f(.11, 0, .11);
  glVertex3f(.11, -1, .11);
  glVertex3f(.11, -1, -.11);
  glVertex3f(.11, 0, -.11);
  glVertex3f(.11, 0, .11);
  glVertex3f(.11, -1, .11);
  glEnd();
}

void zScene() {
  glBegin(GL_LINE_STRIP);
  glColor3f(1.0, 1.0, 1.0);
  glVertex3f(-.11, -1, -.21);
  glVertex3f(-.11, 0, -.21);
  glVertex3f(-.11, 0, .11);
  glVertex3f(-.11, -1, .11);
  glVertex3f(-.11, -1, -.21);
  glVertex3f(.11, -1, -.21);
  glVertex3f(.11, 0, -.21);
  glVertex3f(-.11, 0, -.21);
  glVertex3f(-.11, -1, -.21);
  glEnd();
  glBegin(GL_LINE_STRIP);
  glColor3f(1.0, 1.0, 1.0);
  glVertex3f(.11, -1, .11);
  glVertex3f(-.11, -1, .11);
  glVertex3f(-.11, 0, .11);
  glVertex3f(.11, 0, .11);
  glVertex3f(.11, -1, .11);
  glVertex3f(.11, -1, -.21);
  glVertex3f(.11, 0, -.21);
  glVertex3f(.11, 0, .11);
  glVertex3f(.11, -1, .11);
  glEnd();
}

void noBound() {
  glBegin(GL_LINE_STRIP);
  glColor3f(1.0, 1.0, 1.0);
  glVertex3f(-maxX, -1, -maxZ);
  glVertex3f(-maxX, maxY, -maxZ);
  glVertex3f(-maxX, maxY, maxZ);
  glVertex3f(-maxX, -1, maxZ);
  glVertex3f(-maxX, -1, -maxZ);
  glVertex3f(maxX, -1, -maxZ);
  glVertex3f(maxX, maxY, -maxZ);
  glVertex3f(-maxX, maxY, -maxZ);
  glVertex3f(-maxX, -1, -maxZ);
  glEnd();
  glBegin(GL_LINE_STRIP);
  glColor3f(1.0, 1.0, 1.0);
  glVertex3f(maxX, -1, maxZ);
  glVertex3f(-maxX, -1, maxZ);
  glVertex3f(-maxX, maxY, maxZ);
  glVertex3f(maxX, maxY, maxZ);
  glVertex3f(maxX, -1, maxZ);
  glVertex3f(maxX, -1, -maxZ);
  glVertex3f(maxX, maxY, -maxZ);
  glVertex3f(maxX, maxY, maxZ);
  glVertex3f(maxX, -1, maxZ);
  glEnd();
}

int camera = 1;
int numcameras = 3;
unsigned char* myBMP;
float angle = 0;
/* Display is updated/rendered here. */
void displayFunc() {

  // Clear Color and Depth Buffers
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Reset transformations
  glLoadIdentity();
  // Set the camera
  
  if(camera==1){
  gluLookAt(	0.0f, 0.3f, 2.0f,
		0.0f, -.5f,  0.0f,
		0.0f, 1.0f,  0.0f);
  }
  else if(camera == 2){
  gluLookAt(	0.0f, -.7f, 2.0f,
		0.0f, -.7f,  0.0f,
		0.0f, 1.0f,  0.0f);
  }
  // new settings just to see obj file. can change later.
  else if(camera == 3){
  gluLookAt(	0.0f, 2.3f, 1.5f,
		0.0f, -.5f,  0.0f,
		0.0f, 1.0f,  0.0f);
  }


  glRotatef(angle, 0.0f, 1.0f, 0.0f);
  glColor4f(1.0, 1.0, 1.0, 1.0);
  glBegin(GL_TRIANGLES);
  glVertex3f( -1.0f, -1.0f, -1.0f);
  glVertex3f( 1.0f, -1.0f, -1.0);
  glVertex3f( 1.0f, -1.0f, 1.0);
  glEnd();

  glBegin(GL_TRIANGLES);
  glVertex3f( 1.0f, -1.0f, 1.0f);
  glVertex3f( -1.0f, -1.0f, 1.0);
  glVertex3f( -1.0f, -1.0f, -1.0);
  glEnd();

  angle+=0.3f;

  if (bBoxScene) {
    boxScene();
  }
  else if (bxScene) {
    xScene();
  }
  else if (bzScene) {
    zScene();
  }
  else if (bNoBound) {
    noBound();
  }
  glClearDepth(1);
  /*glEnable(GL_DEPTH_TEST);
  glEnable(GL_LIGHTING);
  glEnable(GL_COLOR_MATERIAL);

  GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat mat_shininess[] = { 50.0 };
  GLfloat light_position[] = { 1.0, 1.0, 1.0, 0.0 };
  glClearColor (0.0, 0.0, 0.0, 0.0);
  glShadeModel (GL_SMOOTH);

  glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
  glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);


  glEnable(GL_LIGHT0);*/

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_LIGHTING);
  glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_LIGHT0);
  float lightPos0[4] = { 0, 0, 5, 1 };
  float dif0[4] = { .6, .6, .6, 1 };
  float amb0[4] = { .1, .1, .1, 1 };
  float spec0[4] = { .6, .6, .6, 1 };
  float shine0[1] = {1};
  glLightfv(GL_LIGHT0, GL_POSITION, lightPos0);
  glLightfv(GL_LIGHT0, GL_AMBIENT, amb0);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, dif0);
  glLightfv(GL_LIGHT0, GL_SPECULAR, spec0);
  glLightfv(GL_LIGHT0, GL_SHININESS, shine0);
  glEnable(GL_LIGHT1);
  float lightPos1[4] = {0, 2, -5, 1};
  glLightfv(GL_LIGHT1, GL_POSITION, lightPos1);


  glBegin(GL_POINTS);
  //glColor4f(0.0, 0.2, 1.0, 1.0);
  Particle* p;
  unsigned int i = 0;
  for (i = 0; i < pSystem.getParticles().size()/2; i++) {
    p = pSystem.getParticle(i);
	glColor4f(p->getColor().getX(),p->getColor().getY(),p->getColor().getZ(), 1.0);
    glVertex3f(p->getPosition().getX(), p->getPosition().getY(), p->getPosition().getZ());
  }
  glEnd();
  glBegin(GL_POINTS);
  //glColor4f(0.2, 0.7, 0.0, 1.0);
  for (; i < pSystem.getParticles().size(); i++) {
    p = pSystem.getParticle(i);
	glColor4f(p->getColor().getX(),p->getColor().getY(),p->getColor().getZ(), 1.0);
    glVertex3f(p->getPosition().getX(), p->getPosition().getY(), p->getPosition().getZ());
  }
  glEnd();
  if (currentTime < totalTime) {
    pSystem.update(dt);
    currentTime += dt;
  } else {
    exit(1);
  }
  //glDisable(GL_LIGHTING);
  glFlush();
  if(record) {
    recorder.RecordFrame();
  }
  glutSwapBuffers();
}

void displayFunc1() {

	  // Clear Color and Depth Buffers
	  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	  // Reset transformations
	  glLoadIdentity();
	  // Set the camera

	  if(camera==1){
	  gluLookAt(	0.0f, 0.3f, 2.0f,
			0.0f, -.5f,  0.0f,
			0.0f, 1.0f,  0.0f);
	  }
	  else if(camera == 2){
	  gluLookAt(	0.0f, -.7f, 2.0f,
			0.0f, -.7f,  0.0f,
			0.0f, 1.0f,  0.0f);
	  }
	  // new settings just to see obj file. can change later.
	  else if(camera == 3){
	  gluLookAt(	0.0f, 2.3f, 1.5f,
			0.0f, -.5f,  0.0f,
			0.0f, 1.0f,  0.0f);
	  }


	  glRotatef(angle, 0.0f, 1.0f, 0.0f);
	  glColor4f(1.0, 1.0, 1.0, 1.0);
	  glBegin(GL_TRIANGLES);
	  glVertex3f( -1.0f, -1.0f, -1.0f);
	  glVertex3f( 1.0f, -1.0f, -1.0);
	  glVertex3f( 1.0f, -1.0f, 1.0);
	  glEnd();

	  glBegin(GL_TRIANGLES);
	  glVertex3f( 1.0f, -1.0f, 1.0f);
	  glVertex3f( -1.0f, -1.0f, 1.0);
	  glVertex3f( -1.0f, -1.0f, -1.0);
	  glEnd();

	  angle+=0.3f;

	  if (bBoxScene) {
	    boxScene();
	  }
	  else if (bxScene) {
	    xScene();
	  }
	  else if (bzScene) {
	    zScene();
	  }
	  else if (bNoBound) {
	    noBound();
	  }
	  glClearDepth(1);
	  /*glEnable(GL_DEPTH_TEST);
	  glEnable(GL_LIGHTING);
	  glEnable(GL_COLOR_MATERIAL);

	  GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
	  GLfloat mat_shininess[] = { 50.0 };
	  GLfloat light_position[] = { 1.0, 1.0, 1.0, 0.0 };
	  glClearColor (0.0, 0.0, 0.0, 0.0);
	  glShadeModel (GL_SMOOTH);

	  glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	  glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
	  glLightfv(GL_LIGHT0, GL_POSITION, light_position);


	  glEnable(GL_LIGHT0);*/

	  glEnable(GL_DEPTH_TEST);
	  glEnable(GL_LIGHTING);
	  glEnable(GL_COLOR_MATERIAL);
	    glEnable(GL_LIGHT0);
	  float lightPos0[4] = { 0, 0, 5, 1 };
	  float dif0[4] = { .6, .6, .6, 1 };
	  float amb0[4] = { .1, .1, .1, 1 };
	  float spec0[4] = { .6, .6, .6, 1 };
	  float shine0[1] = {1};
	  glLightfv(GL_LIGHT0, GL_POSITION, lightPos0);
	  glLightfv(GL_LIGHT0, GL_AMBIENT, amb0);
	  glLightfv(GL_LIGHT0, GL_DIFFUSE, dif0);
	  glLightfv(GL_LIGHT0, GL_SPECULAR, spec0);
	  glLightfv(GL_LIGHT0, GL_SHININESS, shine0);
	  glEnable(GL_LIGHT1);
	  float lightPos1[4] = {0, 2, -5, 1};
	  glLightfv(GL_LIGHT1, GL_POSITION, lightPos1);


	  mesh.updateColors();
	  vector<float> triangles, normals;
	  mesh.marchingCubes(triangles, normals);

	  glBegin(GL_TRIANGLES);
	  glColor4f(0.1, 0.3, 0.8, 0.7);
	  for (int j = 0; j < (int) triangles.size(); j+=3) {
		glNormal3f(normals[j], normals[j+1], normals[j+2]);
		glVertex3f(triangles[j], triangles[j+1], triangles[j+2]);
	  }
	  glEnd();

	  if (currentTime < totalTime) {
	    pSystem.update(dt);
	    currentTime += dt;
	  } else {
	    exit(1);
	  }
	  //glDisable(GL_LIGHTING);
	  glFlush();
	  if(record) {
	    recorder.RecordFrame();
	  }
	  glutSwapBuffers();
}

void keyboard(unsigned char key, int x, int y) {
  Particle p1;
  switch(key) {
    case ' ': // allow spacebar to end the program
      exit(0);
      break;
    case 'b':
      maxX = 1;
      maxY = 0;
      maxZ = 1;
      bBoxScene = false;
      bxScene = false;
      bzScene = false;
      bNoBound = true;
      pSystem.setBoundaries(-maxX, -1, -maxZ, maxX, maxY, maxZ);
      break;
    case 'p':
      p1 = Water();
      p1.setPosition(Point3D(0,.1,0));
      pSystem.addParticle(p1);
      break;
	case 'o':
		for(float i = 0.; i < 5; i++){
			for(float j = 0.; j < 5; j++){
				for(float k = 0.; k < 5; k++){
					p1 = Mucus();
					p1.setPosition(Point3D(i/50,j/50-.2,k/50));
					pSystem.addParticle(p1);
				}
			}
		}
      break;
	  	case 'w':
		for(float i = 0.; i < 5; i++){
			for(float j = 0.; j < 5; j++){
				for(float k = 0.; k < 5; k++){
					p1 = Water();
					p1.setPosition(Point3D(i/50,j/50-.2,k/50));
					pSystem.addParticle(p1);
				}
			}
		}
      break;
    case 'x': // opens up x portion of boundary
      bBoxScene = false;
      bxScene = true;
      bzScene = false;
      bNoBound = false;
      pSystem.setBoundaries(-1, -1, -.1, .1, 0, .1);
      break;
    case 'z': // open up z portion of the boundary
      bBoxScene = false;
      bxScene = false;
      bzScene = true;
      bNoBound = false;
      pSystem.setBoundaries(-.1, -1, -.2, .1, 0, .1);
      break;
	case 'c': // change camera
		camera++;
		if (camera > numcameras)
			camera = 1;
		break;
    default:
      break;
  }
  glutPostRedisplay();
}

void changeSize(int w, int h) {

  // Prevent a divide by zero, when window is too short
  // (you cant make a window of zero width).
  if (h == 0)
    h = 1;

  float ratio =  w * 1.0 / h;

  // Use the Projection Matrix
  glMatrixMode(GL_PROJECTION);

  // Reset Matrix
  glLoadIdentity();

  // Set the viewport to be the entire window
  glViewport(0, 0, w, h);

  // Set the correct perspective.
  gluPerspective(45.0f, ratio, 0.1f, 100.0f);

  // Get Back to the Modelview
  glMatrixMode(GL_MODELVIEW);
}

int main(int argc, char** argv) {
  pSystem = ParticleSystem(Vector(0,-9.8, 0));
  determineFunction(argc, argv);
  createParticles();
  pSystem.initialize(dt);
  // dont need this now that we have opengl.
  // but we can use this if we don't want to display...i guess
  /*for(float j = 0; j < totalTime; j+=dt) {*/
  //cout << "//===========================================//" << endl
  //<< "// totalTime: " << j << "                                 //" << endl
  //<< "//===========================================//" << endl;
  //pSystem.update(dt);
  /*}*/
  glutInit(&argc, argv);
  glutInitWindowSize(frameWidth, frameHeight);
  glutInitWindowPosition(0, 0);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  glutCreateWindow("CS184 Final");
  glClearColor(0.0, 0.0, 0.0, 1.0);
  glEnable(GL_POINT_SMOOTH);
  glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);
  glPointSize(5.0f);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_BLEND);
  glutKeyboardFunc(keyboard);

  if (useCubes) {
      mesh = Mesh(32, 33, 18, 1.0, -1.0, .5, -1.5, .5, -.5);
	  glutDisplayFunc(displayFunc1);
	  glutIdleFunc(displayFunc1);
  } else {
	  glutDisplayFunc(displayFunc);
	  glutIdleFunc(displayFunc);
  }
  glutReshapeFunc(changeSize);
  glutMainLoop();
  return 0;
}


