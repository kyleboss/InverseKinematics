using Eigen::MatrixXd;
using Eigen::Vector3d;

class Segment {
  public:
    Segment(float length) : length(length) {
        transMatrix = Eigen::Matrix<double,3,3>::Identity();
<<<<<<< HEAD
        testMatrix = Eigen::Matrix<double,3,3>::Identity();
=======
>>>>>>> 6124439b3a37ecdd6f81c428f34030ccd51699f8
        numSegments += 1;
        jointLoc = Vector3d(totalLength,0,0);
        // oldLoc = Vector3d(totalLength,0,0);
        end = Vector3d(totalLength+length,0,0);
        // oldEnd = Vector3d(totalLength+length,0,0);
        totalLength += length;
    }
    Segment() {}
    Vector3d jointLoc; //also the start drawing location
    Vector3d end;
    // Vector3d oldEnd;
    // Vector3d oldLoc;
    Eigen::Matrix<double,3,3> transMatrix;
<<<<<<< HEAD
    Eigen::Matrix<double,3,3> testMatrix;
=======
>>>>>>> 6124439b3a37ecdd6f81c428f34030ccd51699f8
    float length;
    static int numSegments;
    static int totalLength;
};
int Segment::numSegments = 0;
int Segment::totalLength = 0;