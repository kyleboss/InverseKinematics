using Eigen::MatrixXd;
using Eigen::Vector3d;

class Segment {
  public:
    Segment(float length) : length(length) {
        transMatrix = Eigen::Matrix<double,3,3>::Identity();
        oldTransMatrix = Eigen::Matrix<double,3,3>::Identity();
        numSegments += 1;
        jointLoc = Vector3d(totalLength,0,0);
        oldLoc = Vector3d(totalLength,0,0);
        end = Vector3d(totalLength+length,0,0);
        oldEnd = Vector3d(totalLength+length,0,0);
        totalLength += length;
    }
    Segment() {}
    Vector3d jointLoc; //also the start drawing location
    Vector3d end;
    Vector3d oldEnd;
    Vector3d oldLoc;
    Eigen::Matrix<double,3,3> transMatrix;
    Eigen::Matrix<double,3,3> oldTransMatrix;
    float length;
    static int numSegments;
    static int totalLength;
};
int Segment::numSegments = 0;
int Segment::totalLength = 0;