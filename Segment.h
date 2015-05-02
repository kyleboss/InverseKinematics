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
    float length;
    static int numSegments;
    static std::vector<Segment> * segments;
};
int Segment::numSegments = 0;
