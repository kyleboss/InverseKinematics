#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <sys/time.h>
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#include <time.h>
#include <math.h>
#include "Eigen/Dense"
#include "Segment.h"

using Eigen::MatrixXd;
using Eigen::Vector3d;
using Eigen::Affine3d;
using Eigen::AngleAxisd;
using Eigen::Translation3d;
using namespace std;
Segment * youngestSeg;
Segment * rootSeg;
float acceptableDistance;
Vector3d goal;
std::vector<Segment> segments = std::vector<Segment>();

//*********************************************************
// distanceBetween
// Returns magnitude of two vectors.
//*********************************************************
float distanceBetween(Vector3d point1, Vector3d point2) {
  return sqrt(pow(point1[0]-point2[0],2)+pow(point1[0]-point2[0],2)+pow(point1[0]-point2[0],2));
}

//*********************************************************
// getEndPoint
// Determines the end-point location of the segment that
// is index-amount of hops away from the root. If no index
// is provided, it will return the end-point of the segment
// farthest away from the root. 
// http://stackoverflow.com/questions/10115354/inverse-kinematics-with-opengl-eigen3-unstable-jacobian-pseudoinverse
//*********************************************************
Vector3d getEndPoint(int index = Segment::numSegments, bool draw = false) {
  Vector3d prevEndPoint;
  prevEndPoint = Vector3d(0,0,0);
  if (draw) {
    glColor3f(1,1,1);
    glBegin(GL_POINTS);
  }
  Vector3d endPoint = Vector3d(0,0,0);
  for (int i = 0; i<index && i<Segment::numSegments; i++) {
    Segment currentSegment    = segments[i];
    Vector3d rad              = M_PI*currentSegment.rot/180;
    AngleAxisd xRot           = AngleAxisd(rad[0], Vector3d(-1, 0, 0));
    AngleAxisd yRot           = AngleAxisd(rad[1], Vector3d(0, -1, 0));
    AngleAxisd zRot           = AngleAxisd(rad[2], Vector3d(0, 0, -1));
    Translation3d translation = Translation3d(Vector3d(0, 0, currentSegment.length));
    endPoint                  = ((Affine3d) xRot*yRot*zRot*translation)*endPoint;
    if (draw) {
      glVertex3f(prevEndPoint[0], prevEndPoint[1], prevEndPoint[2]);
      glVertex3f(endPoint[0],endPoint[1],endPoint[2]);
      cout << "for segment of len " << currentSegment.length << " prev is " 
           << prevEndPoint << " and curr is " << endPoint << endl;
    }
    if (draw) {
      prevEndPoint = endPoint;
    }
  }
  if (draw) {
    glEnd();
  }
  return endPoint;
}

//*********************************************************
// computeJacobian
// Computes the jacobian matrix from the current segments &
// endpoints.
// http://www.dandiggins.co.uk/iksolver-3.html
//*********************************************************
MatrixXd computeJacobian() {
  MatrixXd jacobian = MatrixXd(3,3*Segment::numSegments);
  Vector3d xVec     = Vector3d(1,0,0);
  Vector3d yVec     = Vector3d(0,1,0);
  Vector3d zVec     = Vector3d(0,0,1);
  for (int i=0; i<Segment::numSegments; i++) {
    Vector3d endPoint   = getEndPoint(i);
    Vector3d difference = goal-endPoint;
    Vector3d xCol       = xVec.cross(difference);
    Vector3d yCol       = yVec.cross(difference);
    Vector3d zCol       = zVec.cross(difference);
    jacobian.col(3*i+0) = xCol;
    jacobian.col(3*i+1) = yCol;
    jacobian.col(3*i+2) = zCol;
  }
  return jacobian;
}

//*********************************************************
// computePseudoInverse
// Computes the pseudoinverse of a given matrix.
//*********************************************************
MatrixXd computePseudoInverse(MatrixXd originalMatrix) {
  return originalMatrix.transpose()*((originalMatrix*originalMatrix.transpose()).inverse());
}

//********************************************************
// inverseKinematicsSolver
// Solves the Inverse Kinematics Problem.
//*********************************************************
void inverseKinematicsSolver() {
  Vector3d endPoint     = getEndPoint();
  float distanceToGoal  = distanceBetween(endPoint, goal);
  double lambda         = 0.1;
  while (distanceToGoal > acceptableDistance) {
    MatrixXd jacobian       = computeJacobian();
    MatrixXd pseudoJacobian = computePseudoInverse(jacobian);
    distanceToGoal          = distanceBetween(endPoint, goal);
    Vector3d addToRots      = pseudoJacobian*lambda*distanceToGoal;
    for (int i = 0; i<Segment::numSegments; i++) {
      segments[i].rot += addToRots;
    } 
    float newDistanceToGoal = distanceBetween(endPoint, goal);
    if (distanceToGoal < newDistanceToGoal) {
      lambda*=.5;
    }
  }
}


//****************************************************
// Some Classes
//****************************************************

class Viewport;

class Viewport {
  public:
    int w, h; // width and height
};


//****************************************************
// Global Variables
//****************************************************
Viewport  viewport;




//****************************************************
// Simple init function
//****************************************************
void initScene(){
  GLfloat black[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat yellow[] = { 1.0, 1.0, 0.0, 1.0 };
  GLfloat cyan[] = { 0.0, 1.0, 1.0, 1.0 };
  GLfloat white[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat direction[] = { 1.0, 1.0, -1.0, 0.0 };

  glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, cyan);
  glMaterialfv(GL_FRONT, GL_SPECULAR, white);
  glMaterialf(GL_FRONT, GL_SHININESS, 30);

  glLightfv(GL_LIGHT0, GL_AMBIENT, black);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, yellow);
  glLightfv(GL_LIGHT0, GL_SPECULAR, white);
  glLightfv(GL_LIGHT0, GL_POSITION, direction);

  glEnable(GL_LIGHTING);                
  glEnable(GL_LIGHT0);  
  glEnable(GL_DEPTH_TEST);
}


//****************************************************
// reshape viewport if the window is resized
//****************************************************
void myReshape(int w, int h) {
  viewport.w = w;
  viewport.h = h;

  glViewport (0,0,viewport.w,viewport.h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

}

//****************************************************
// function that does the actual drawing of stuff
//***************************************************
void myDisplay() {

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-3.5, 3.5, -3.5, 3.5, 5, -5);
 
  // Start drawing
  getEndPoint(Segment::numSegments, true);

  glFlush();
  glutSwapBuffers();          // swap buffers (we earlier set double buffer)
}

void handle(unsigned char key, int x, int y) {
  switch (key) {
    case 32: //space
      exit(0);
      break;
  }
  glutPostRedisplay();
}



int main (int argc, char *argv[]) {

  MatrixXd m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);
  std::cout << m << std::endl;

  Segment a = Segment(1);
  Segment b = Segment(4);
  Segment c = Segment(2);
  segments.push_back(a);
  segments.push_back(b);
  segments.push_back(c);


  //This initializes glut
  glutInit(&argc, argv);


  //This tells glut to use a double-buffered window with red, green, and blue channels 
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

  // Initalize theviewport size
  viewport.w = 700;
  viewport.h = 700;

  //The size and position of the window
  glutInitWindowSize(viewport.w, viewport.h);
  glutInitWindowPosition(0,0);
  glutCreateWindow(argv[0]);

  initScene();              // quick function to set up scene
  glutDisplayFunc(myDisplay);       // function to run when its time to draw something
  glutReshapeFunc(myReshape);       // function to run when the window gets resized
  glutKeyboardFunc(handle); //exit on space
  glutMainLoop();             // infinite loop that will keep drawing and resizing
  

  // and whatever else
  return 0;
}