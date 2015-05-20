#include <iostream>
#include <vector>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <sphep/generator.hh>
#include <Eigen/Dense>

template <class V>
void writeRaw(const std::string& name, const std::vector<V>& points)
{
  std::ofstream stream(name);
  for (auto& p : points) {
    for (unsigned int i = 0; i < p.rows(); ++i) {
      stream << (i > 0 ? " " : "") << p(i);
    }
    stream << "\n";
  }
}

int main(int argc, char** argv)
{
  if (argc != 2) {
    std::cerr << "please provide the number of points N. call:\n" << argv[0] << " <N>\n";
    return -1;
  }
  typedef double T;
  const int dim = 3;
  typedef Eigen::Matrix<T, dim, 1> Vector;
  std::vector<Vector> points;
  unsigned int N = std::stoi(argv[1]);
  sphep::generatePoints(N, points);
  writeRaw("points.txt", points);
}
