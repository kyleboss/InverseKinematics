using Eigen::MatrixXd;
using Eigen::Vector3d;

class Segment {
  public:
    Segment(float length) : length(length) {
        transMatrix = Eigen::Matrix<double,3,3>::Identity();
        numSegments += 1;
        jointLoc = Vector3d(totalLength,0,0);
        end = Vector3d(totalLength+length,0,0);
        old_end = end;
        old_jointLoc = jointLoc;
        totalLength += length;
    }
    Segment() {}
    Vector3d jointLoc; //also the start drawing location
    Vector3d end;
    Vector3d old_jointLoc;
    Vector3d old_end;
    Eigen::Matrix<double,3,3> transMatrix;
    float length;
    static int numSegments;
    static int totalLength;
};
int Segment::numSegments = 0;
int Segment::totalLength = 0;