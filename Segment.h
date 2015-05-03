using Eigen::MatrixXd;
using Eigen::Vector3d;

class Segment {
  public:
    Segment(float length) : length(length) {
        rot = Vector3d(0,0,0);
        numSegments += 1;
    }
    Segment() {}
    Vector3d rot;
    Vector3d start;
    Vector3d end;
    float length;
    static int numSegments;
};
int Segment::numSegments = 0;
