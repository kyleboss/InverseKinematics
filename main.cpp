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
using Eigen::VectorXd;
using Eigen::Affine3d;
using Eigen::AngleAxisd;
using Eigen::Translation3d;
using namespace std;
Segment * youngestSeg;
Segment * rootSeg;
float acceptableDistance      = .001;
Vector3d goal                 = Vector3d(2, 2, 0);
std::vector<Segment> segments = std::vector<Segment>();

//*********************************************************
// distanceBetween
// Returns magnitude of two vectors.
//*********************************************************
float distanceBetween(Vector3d point1, Vector3d point2) {
  return sqrt(pow(point1[0]-point2[0],2)+pow(point1[1]-point2[1],2)+pow(point1[2]-point2[2],2));
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
  Vector3d prevEndPoint = Vector3d(0,0,0);
  if (draw) {
    glPointSize(6);
    glLineWidth(6);
    glBegin(GL_LINES);
  }
  Vector3d endPoint = Vector3d(0,0,0);
  for (int i = 0; i<index && i<Segment::numSegments; i++) {
    Segment currentSegment    = segments[i];
    Vector3d rad              = M_PI*currentSegment.rot/180;
    AngleAxisd xRot           = AngleAxisd(rad[0], Vector3d(-1, 0, 0));
    AngleAxisd yRot           = AngleAxisd(rad[1], Vector3d(0, -1, 0));
    AngleAxisd zRot           = AngleAxisd(rad[2], Vector3d(0, 0, -1));
    Translation3d translation = Translation3d(Vector3d(currentSegment.length, 0, 0));
    endPoint                  = ((Affine3d) xRot*yRot*zRot*translation)*endPoint;
    if (draw) {
      if (i==0) {
        glColor3f(1,0,0);
      } else if(i==1) {
        glColor3f(0,1,0);
      } else if(i==2) {
        glColor3f(0,0,1);
      } else if(i==3) {
        glColor3f(1,0,1);
      }
      cout << "distance: " << distanceBetween(prevEndPoint, endPoint) << endl;  
      glVertex3d(prevEndPoint[0], prevEndPoint[1], prevEndPoint[2]);
      glVertex3d(endPoint[0], endPoint[1], endPoint[2]);
    }
    if (draw) {
      prevEndPoint = endPoint;
    }
  }
  if (draw) {
    glEnd();
    glColor3f(1,1,1);
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
    Vector3d endPoint   = getEndPoint(i+1);
    Vector3d difference = goal-endPoint;
    Vector3d xCol = Vector3d(0,0,0);
    Vector3d yCol = Vector3d(0,0,0);
    Vector3d zCol = Vector3d(0,0,0);
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
  float newDistanceToGoal;
  MatrixXd jacobian;
  MatrixXd pseudoJacobian;
  VectorXd addToRots;
  int numCalcsBeforeGiveUp = 10000;
  int numCalcs = 0;
  cout << "Distance: " << distanceToGoal << endl;
  cout << "endPoint: " << endPoint << endl;
  cout << "goal: " << goal << endl;
  while (distanceToGoal > acceptableDistance && numCalcs < numCalcsBeforeGiveUp) {
    numCalcs++;
    cout << "Distance: " << distanceToGoal << endl;
    jacobian       = computeJacobian();
    pseudoJacobian = computePseudoInverse(jacobian);
    distanceToGoal          = distanceBetween(endPoint, goal);
    addToRots      = pseudoJacobian*lambda*(goal-endPoint);
    for (int i = 0; i<Segment::numSegments; i++) {
      segments[i].rot[0] = fmod((segments[i].rot[0] + addToRots[i*3+0]), 360);
      segments[i].rot[1] = fmod((segments[i].rot[1] + addToRots[i*3+1]), 360);
      segments[i].rot[2] = fmod((segments[i].rot[2] + addToRots[i*3+2]), 360);
      // cout << "segments[i[.rot: " << jacobian << endl;
    } 
    endPoint     = getEndPoint();
    newDistanceToGoal = distanceBetween(endPoint, goal);
    if (distanceToGoal < newDistanceToGoal) {
      lambda*=.5;
    }
  }
}

#define PI 3.14159265

using namespace std;

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
  viewport.w = w;
  viewport.h = h;

  glViewport(0,0,viewport.w,viewport.h);// sets the rectangle that will be the window
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();                // loading the identity matrix for the screen

  //----------- setting the projection -------------------------
  // glOrtho sets left, right, bottom, top, zNear, zFar of the chord system


  // glOrtho(-1, 1 + (w-400)/200.0 , -1 -(h-400)/200.0, 1, 1, -1); // resize type = add
  // glOrtho(-w/400.0, w/400.0, -h/400.0, h/400.0, 1, -1); // resize type = center

  glOrtho(-9, 9, -9, 9, 9, -9);    // resize type = stretch

  //------------------------------------------------------------
}


//****************************************************
// sets the window up
//****************************************************
void initScene(){
  glClearColor(0.0f, 0.0f, 0.0f, 0.0f); // Clear to black, fully transparent

  myReshape(viewport.w,viewport.h);
}

void makeCircle(float rad, float horizOffset, float vertOffset, float distortY) {
  const  float circ = 2*PI;
  for (float angle = 0; angle < circ; angle += 0.01) {
    glVertex3f(rad*cos(angle)+horizOffset, distortY*rad*sin(angle)+vertOffset, 0.0f);
  }
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

  glClear(GL_COLOR_BUFFER_BIT);                // clear the color buffer (sets everything to black)

  glMatrixMode(GL_MODELVIEW);                  // indicate we are specifying camera transformations
  glLoadIdentity();                            // make sure transformation is "zero'd"

  //----------------------- code to draw objects --------------------------

  glColor3f(0.75f,1.0f,0.0f);
  Segment a = Segment(1);
  Segment b = Segment(1);
  // Segment c = Segment(1);
  // Segment d = Segment(1);
  segments.push_back(a);
  segments.push_back(b);
  // segments.push_back(c);
  // segments.push_back(d);
  inverseKinematicsSolver();
  getEndPoint(Segment::numSegments, true);
  glBegin(GL_POINTS);
  glPointSize(10);
  glVertex3f(goal[0], goal[1], goal[2]);
  glVertex3f(0, 0, 0);
  glEnd();


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


//****************************************************
// the usual stuff, nothing exciting here
//****************************************************
int main(int argc, char *argv[]) {
  //This initializes glut
  glutInit(&argc, argv);

  //This tells glut to use a double-buffered window with red, green, and blue channels 
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);

  // Initalize theviewport size
  viewport.w = 400;
  viewport.h = 400;

  //The size and position of the window
  glutInitWindowSize(viewport.w, viewport.h);
  glutInitWindowPosition(0, 0);
  glutCreateWindow("CS184!");

  initScene();              // quick function to set up scene
  glutDisplayFunc(myDisplay);       // function to run when its time to draw something
  glutReshapeFunc(myReshape);       // function to run when the window gets resized
  glutKeyboardFunc(handle); //exit on space
  glutIdleFunc(myFrameMove);                   // function to run when not handling any other task

  glutMainLoop();             // infinite loop that will keep drawing and resizing
 
  return 0;
}