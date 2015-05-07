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
#include "Eigen/SVD"
#include "Segment.h"

#define PI 3.14159265

using Eigen::JacobiSVD;
using Eigen::ComputeThinU;
using Eigen::ComputeThinV;
using Eigen::MatrixXd;
using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::Affine3d;
using Eigen::AngleAxisd;
using Eigen::Translation3d;
using namespace std;

Segment * youngestSeg;
Segment * rootSeg;
int timeCount = 0;
float acceptableDistance      = .1;
Vector3d goal                 = Vector3d(2, 0, 0);
Vector3d realGoal                 = Vector3d(2, 0, 0);

std::vector<Segment *> segments = std::vector<Segment *>();

//*********************************************************
// distanceBetween
// Returns magnitude of two vectors.
//*********************************************************
float distanceBetween(Vector3d point1, Vector3d point2) {
  return sqrt(pow(point1[0]-point2[0],2)+pow(point1[1]-point2[1],2)+pow(point1[2]-point2[2],2));
}

//*********************************************************
// changeColor
// Changes the color of the openGL outputs.
//*********************************************************
void changeColor(float r, float g, float b) {
  glColor3f(r,g,b);
}

//*********************************************************
// alterColorForDebugging
// Helper function to better understand the IK Solver from
// within the getEndPoint function.
//*********************************************************
void alterColorForDebugging(int i, Vector3d prevEndPoint, Vector3d endPoint) {
  glPointSize(6);
  glLineWidth(6);
  glBegin(GL_LINES);
  if (i==0) glColor3d(1,0,0);
  if (i==1) glColor3d(0,1,0);
  if (i==2) glColor3d(0,0,1);
  if (i==3) glColor3d(1,0,1);
  glVertex3d(prevEndPoint[0], prevEndPoint[1], prevEndPoint[2]);
  glVertex3d(endPoint[0], endPoint[1], endPoint[2]);
  glEnd();
}

//*********************************************************
// getEndPoint
// Determines the end-point location of the segment that
// is index-amount of hops away from the root. If no index
// is provided, it will return the end-point of the segment
// farthest away from the root. 
//*********************************************************
Vector3d getEndPoint(int index = Segment::numSegments, bool draw = false, bool test = false) {
  Vector3d endPoint = Vector3d(0,0,0);
  Vector3d prevEndPoint = Vector3d(0,0,0);

  for (int i = 0; i<index && i<Segment::numSegments; i++) {
    if (!test) {
      endPoint += segments[i]->transMatrix*Vector3d(segments[i]->length,0,0);
      segments[i]->end = endPoint;
      segments[i]->jointLoc = prevEndPoint;
    } else {
      endPoint += segments[i]->testMatrix*Vector3d(segments[i]->length,0,0);
    }
    if (draw) alterColorForDebugging(i, prevEndPoint, endPoint);
    prevEndPoint = endPoint;
  }
  return endPoint;
}

//*********************************************************
// computeJacobian
// Computes the jacobian matrix from the current segments &
// endpoints.
//*********************************************************
MatrixXd computeJacobian() {
  MatrixXd jacobian = MatrixXd(3,3*Segment::numSegments);
  Vector3d xVec     = Vector3d(1,0,0);
  Vector3d yVec     = Vector3d(0,1,0);
  Vector3d zVec     = Vector3d(0,0,1);
  Vector3d xCol, yCol, zCol, joint, difference;
  Vector3d endEffector = segments[Segment::numSegments-1]->end;
  //need to get into world coordinate space
  for (int i=0; i<Segment::numSegments; i++) {
    Segment * currentSegment = segments[i];
    joint = currentSegment->jointLoc;
    difference  = endEffector-joint; //difference = end effector - curr joint
    xCol        = Vector3d(0,0,0);
    yCol        = Vector3d(0,0,0);
    zCol        = Vector3d(0,0,0);

    xCol       += xVec.cross(difference);
    yCol       += yVec.cross(difference);
    zCol       += zVec.cross(difference);
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
MatrixXd computePseudoInverse(MatrixXd originalMatrix, Vector3d goal, Vector3d endPoint) {
    JacobiSVD<MatrixXd> svd(originalMatrix,ComputeThinU | ComputeThinV);
    MatrixXd sigma_inverse = MatrixXd::Zero(3,3);
    Vector3d sigma = svd.singularValues();
    for (int i = 0; i < sigma.size(); i++) { //manually insert the values
      if (sigma(i) > 1e-6) sigma_inverse(i,i) = 1/sigma(i);
    }
    MatrixXd inversejacobian = svd.matrixV() * sigma_inverse * svd.matrixU().transpose();
    
    return inversejacobian;
}

//*********************************************************
// updateSegmentRotations
// Updates a segment's degree of rotation given a vector
// of rotational degree values. addToRots - 1x3n
// rotations are in degreeees
//*********************************************************
void updateSegmentRotations(VectorXd addToRots, bool test=false) {
  AngleAxisd rotx;
  AngleAxisd roty;
  AngleAxisd rotz;
  AngleAxisd rot;
  Segment * currentSegment;
  for (int i = 0; i<Segment::numSegments; i++) { //x, y, z
    currentSegment  = segments[i];
    rot = AngleAxisd(addToRots[3*i+0], currentSegment->transMatrix*Vector3d(1,0,0));
    rot = AngleAxisd(addToRots[3*i+1], currentSegment->transMatrix*Vector3d(0,1,0))*rot;
    rot = AngleAxisd(addToRots[3*i+2], currentSegment->transMatrix*Vector3d(0,0,1))*rot;
    if (!test) {
      currentSegment->transMatrix = rot*currentSegment->transMatrix;
    } else {
      currentSegment->testMatrix = rot*currentSegment->transMatrix;
    }
  } 
}

//********************************************************
// inverseKinematicsSolver
// Solves the Inverse Kinematics Problem.
//*********************************************************
void inverseKinematicsSolver() {
  Vector3d endPoint         = getEndPoint();
  float distanceToGoal      = distanceBetween(endPoint, goal);
  double lambda             = 1;
  int numCalcs              = 0;
  float newDistanceToGoal;
  MatrixXd jacobian;
  MatrixXd pseudoJacobian;
  VectorXd addToRots = Vector3d(0,0,0);

  // if (Segment::totalLength < distanceBetween(Vector3d(0,0,0), goal)) { 
  //   goal = goal.normalized() * Segment::totalLength; 
  // }


  while (distanceToGoal > acceptableDistance && numCalcs < 100*Segment::numSegments) {
    numCalcs++;
    jacobian       = computeJacobian();
    pseudoJacobian = computePseudoInverse(jacobian, goal, endPoint);
    addToRots      = pseudoJacobian*(goal - endPoint);
    updateSegmentRotations(addToRots*lambda, true);
    endPoint = getEndPoint(Segment::numSegments, false, true); //correct reupdating?
    newDistanceToGoal = distanceBetween(endPoint, goal);
    // cout << "newDistanceToGoal: " << newDistanceToGoal << endl;
    if (distanceToGoal > newDistanceToGoal) {
      lambda = 1;
      updateSegmentRotations(addToRots*lambda);
      endPoint = getEndPoint();
      newDistanceToGoal = distanceBetween(endPoint, goal);
    } else {
      lambda *= .5;
    }
    endPoint = getEndPoint(Segment::numSegments, true);
    distanceToGoal = distanceBetween(endPoint, goal);
    // glLoadIdentity();
    // glBegin(GL_LINES); 
    // for (int i=0; i<Segment::numSegments; i++) {
    //   glVertex3d(segments[i]->jointLoc[0], segments[i]->jointLoc[1], segments[i]->jointLoc[2]);
    //   glVertex3d(segments[i]->end[0], segments[i]->end[1], segments[i]->end[2]);
    // }
    // glEnd();
  }
}


//****************************************************
// Some Classes
//****************************************************
class Viewport {
  public:
    int w, h; // width and height
};


//****************************************************
// Global Variables
//****************************************************
Viewport    viewport;

//****************************************************
// reshape viewport if the window is resized
//****************************************************
void myReshape(int w, int h) {
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(40.0, GLfloat(w) / GLfloat(h), 1.0, 150.0);
  glMatrixMode(GL_MODELVIEW);
}

//****************************************************
// sets the window up
//****************************************************
void initScene(){
  GLfloat WHITE[] = {1.0,1.0,1.0};
  glEnable(GL_DEPTH_TEST);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, WHITE);
  glLightfv(GL_LIGHT0, GL_SPECULAR, WHITE);
  glMaterialfv(GL_FRONT, GL_SPECULAR, WHITE);
  glMaterialf(GL_FRONT, GL_SHININESS, 30);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glClearColor(0.0f, 0.0f, 0.0f, 0.0f); // Clear to black, fully transparent
  changeColor(0.75f,1.0f,0.0f);
  Segment * a = new Segment(1);
  Segment * b = new Segment(1);
  // Segment * c = new Segment(3);
  // Segment * d = new Segment(4);
  segments.push_back(a);
  segments.push_back(b);
  // segments.push_back(c);
  // segments.push_back(d);
  myReshape(viewport.w,viewport.h);
}

void handle(unsigned char key, int x, int y) {
  switch (key) {
    case 32: //space
      exit(0);
      break;
  }
  glutPostRedisplay();
}

//***************************************************
// function that does the actual drawing
//***************************************************
void myDisplay() {


  //----------------------- code to draw objects --------------------------
  // https://www.opengl.org/discussion_boards/showthread.php/164180-Draw-a-checker-floor
  GLfloat red[] = {1.0,0,0};
  GLfloat orange[] = {1.0,0.4,0};
  GLfloat yellow[] = {1.0,1.0,0};
  GLfloat green[] = {0,1.0,0};
  GLfloat blue[] = {0,0,1.0};
  GLfloat indigo[] = {0.3,0,0.5};
  GLfloat violet[] = {1.0,0,1.0};
  GLfloat white[] = {1.0,1.0,1.0};
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  gluLookAt(5, 2, 8,
            0, 0, 0,
            0.0, 1.0, 0.0);
  GLfloat lightPosition[] = {0, 2, 0, 1};
  glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
  glBegin(GL_QUADS);
  glNormal3d(0, 1, 0);
  for (int x = -4; x < 4 - 1; x++) {
    for (int z = -4; z < 4 - 1; z++) {
      if (abs(x+z)%2 == 0) glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, green);
      if (abs(x+z)%2 == 1) glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, violet);
      glVertex3d(x, 0, z);
      glVertex3d(x+1, 0, z);
      glVertex3d(x+1, 0, z+1);
      glVertex3d(x, 0, z+1);
    }
  }
  glEnd();
  glPointSize(16);
  glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, white);
  glBegin(GL_POINTS);
  // glVertex3d(0,0,0);
  glVertex3f(realGoal[0], realGoal[1], realGoal[2]);
  glEnd();
  inverseKinematicsSolver();
  getEndPoint(Segment::numSegments, true);

  //-----------------------------------------------------------------------

  glFlush();
  glutSwapBuffers();                           // swap buffers (we earlier set double buffer)
}

//****************************************************
// called by glut when there are no messages to handle
//****************************************************
void myFrameMove() {
//   //nothing here for now
// #ifdef _WIN32
//   Sleep(10);                                   //give ~10ms back to OS (so as not to waste the CPU)
// #endif

//   glutPostRedisplay(); // forces glut to call the display function (myDisplay())
}


void timer(int v) {
  glLoadIdentity();
  timeCount++;
  // goal[1] = sin(timeCount)+0;
  // realGoal[1] = sin(timeCount)+0;
  goal[1] = .5*(cos(timeCount))+1;
  realGoal[1] = .5*(cos(timeCount))+1;
  // goal[2] = sin(timeCount)+0;
  // realGoal[1] = .5*(cos(timeCount))+1;
  // goal[0] = 0;
  // goal[1] = 1;
  // goal[2] = 0;
  glutPostRedisplay();
  glutTimerFunc(120, timer, v);
}


//****************************************************
// the usual stuff, nothing exciting here
//****************************************************
int main(int argc, char *argv[]) {
  //This initializes glut
  glutInit(&argc, argv);

  //This tells glut to use a double-buffered window with red, green, and blue channels 
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  viewport.w = 800;
  viewport.h = 800;

  //The size and position of the window
  glutInitWindowSize(viewport.w, viewport.h);
  glutInitWindowPosition(0, 0);
  glutCreateWindow("CS184!");
  glutTimerFunc(100, timer, 0);
  initScene();                                 // quick function to set up scene
  glutKeyboardFunc(handle); //exit on space

  glutDisplayFunc(myDisplay);                  // function to run when its time to draw something
  glutReshapeFunc(myReshape);                  // function to run when the window gets resized
  glutMainLoop();                              // infinite loop that will keep drawing and resizing and whatever else

  return 0;
}