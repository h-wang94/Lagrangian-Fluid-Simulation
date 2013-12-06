/*
 * window.cpp
 *
 *  Created on: Dec 2, 2013
 *      Author: Owner
 */


#include <windows.h>
#include <GL/glut.h>

/* Display is updated/rendered here. */
void displayFunc() {
	glClear(GL_COLOR_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// rendering stuff goes here
	//for each particle
	glBegin(GL_POINTS);
	glColor4f(0.8, 0.2, 0.0, 0.5);
	for (float i = 0.0; i < 1.0; i+=0.1) {
	    glVertex3f(i, 0.0, -1.0);
	}
	glEnd();
	glFlush();
	glutSwapBuffers();
}

int main(int argc, char** argv) {
	glutInit(&argc, argv);
	glutInitWindowSize(600, 600);
	glutInitWindowPosition(5, 5);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	glutCreateWindow("CS184 Final");
	glClearColor(0.0, 0.0, 0.0, 1.0);
	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);
	glPointSize(15.0f);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);

	glutDisplayFunc(displayFunc);
	//glutIdleFunc(displayFunc); if we do user interaction
	glutMainLoop();
	return 0;
}


