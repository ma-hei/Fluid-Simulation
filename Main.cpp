
#include "Grid2.h"
#include <GL/glut.h>
#include <iostream>
#include "math.h"
#define window_width  500
#define window_height 500

using namespace std;

Grid2* grid = new Grid2(60);

void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glColor3f(1,0.5,0.5);
    float velo;

    glBegin(GL_POINTS);
       for (int i=0;i<grid->gridSize*grid->gridSize*4;i++)
       {
    	velo = sqrt(grid->getMarker(i).getVelocity()[0]*grid->getMarker(i).getVelocity()[0]+grid->getMarker(i).getVelocity()[1]*grid->getMarker(i).getVelocity()[1]);
    	glColor3f(velo,0.5,1-velo);
       	glVertex3f(grid->getMarker(i).getPosition()[0]-grid->gridSize/2,grid->getMarker(i).getPosition()[1]-grid->gridSize/2,0);
       }
    glEnd();

    glutSwapBuffers();
}

void reshape(int w, int h) {
	   glViewport (0, 0, (GLsizei) w, (GLsizei) h);
	   glMatrixMode(GL_PROJECTION);
	   glLoadIdentity();
	   glOrtho(-40, 40, -40, 40, -30.0, 30.0);
	   gluLookAt(0,0,0.001,0,0,0,0,10,0);
	   glMatrixMode(GL_MODELVIEW);
	   glLoadIdentity();
}

int counter=0;
void idleFunc(){

		counter++;
		if (counter<100){
		grid->setv(30,4,15);
		}
		grid->advectAll();
		grid->solvePressure();
		grid->updateVelocities();
		grid->moveMarkers();
		//grid->printAllAdvecteduValues();
		//grid->printAlluValues();

	glutPostRedisplay();

}

void test(){

	int counter=0;

	for (int i=0;i<100;i++){
		counter++;
		cout<<counter<<"\n";
		if (counter<200){
		grid->setu(4,5,5);
		}
		grid->advectAll();
		grid->solvePressure();
		grid->updateVelocities();
		grid->moveMarkers();
		grid->printAllAdvecteduValues();
		grid->printAlluValues();
	//glutPostRedisplay();
	}

}

int main(int argc, char** argv) {
	glutInit(&argc, argv);
	glutInitWindowSize(window_width, window_height);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutCreateWindow("water Experiment");
	glutDisplayFunc (display);
	glutReshapeFunc (reshape);
	glutIdleFunc (idleFunc);
	glutMainLoop();
	return 0;

}
