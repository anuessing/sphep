#ifndef SPHEP_GENERATOR_HH
#define SPHEP_GENERATOR_HH

#include <Eigen/Dense>

namespace sphep
{
  // compute force from b to a
  template <class V>
  V computeForce(const V& a, const V& b)
  {
    // compute direction of force a-b
    auto diff = a - b;
    // scale force by norm(diff)^3
    typename V::Scalar diffsq = diff.squaredNorm();
    return diff / (diffsq * std::sqrt(diffsq));
  }

  // update the veloctiy in an explicit euler fashion
  template <class V>
  void updateVelocity(const std::vector<V>& points, typename V::Scalar dt, std::vector<V>& velocity)
  {
    assert(points.size() == velocity.size());

    for (std::size_t i = 0; i < points.size(); ++i) {
      for (std::size_t j = i + 1; j < points.size(); ++j) {
        // compute force from j to i
        V diff = computeForce(points[i], points[j]);
        // update velocity of i and j with dt*diff (using action = reaction)
        velocity[i] += dt * diff;
        velocity[j] -= dt * diff;
      }
    }
  }

  // redistribute the given points on the unit sphere, such that they are
  // more evenly distributed
  template <class V>
  void optimizePointsOnUnitSphere(std::vector<V>& points, unsigned int iterations,
                                  typename V::Scalar dist, typename V::Scalar dt)
  {
    const std::size_t N = points.size();
    std::vector<V> velocity(points.size(), V::Zero());
    for (unsigned int iter = 0; iter < iterations; ++iter) {
      // update velocity (explicit euler)
      updateVelocity(points, dt, velocity);
      // update points (explicit euler)
      typename V::Scalar globalDist = 0.0;
      for (std::size_t i = 0; i < N; ++i) {
        V oldPoint(points[i]);
        // update point
        points[i] += dt * velocity[i];
        // fetch p_i back onto sphere
        points[i] /= points[i].norm();
        // compute distance between old and newpoint
        globalDist += (oldPoint - points[i]).norm();
      }
      // stop iteration if global update is less than threshold
      if (globalDist < dist) {
        break;
      }
    }
  }

  // non-specialized struct for generating random points on the unit sphere
  template <class Vector>
  struct PointOnUnitSphereGenerator;

  // three dimensional spcialization of the random points on unit sphere
  // generator. Uses sphere coordinates.
  template <class T>
  struct PointOnUnitSphereGenerator<Eigen::Matrix<T, 3, 1> >
  {
    typedef Eigen::Matrix<T, 3, 1> Vector;

    Vector operator()() const
    {
      // generate random sphere coordinates between [0,2*PI]
      T theta = 2.0 * M_PI * static_cast<T>(std::rand()) / RAND_MAX;
      T phi = 2.0 * M_PI * static_cast<T>(std::rand()) / RAND_MAX;
      // transform sphere coordinates to cartesian
      Vector p;
      p << std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta);
      return p;
    }
  };

  // perform the operation axpy on a vector
  template <class Vector>
  struct ScaleAndTranslateFunctor
  {
    typedef typename Vector::Scalar Scalar;

    ScaleAndTranslateFunctor(const Scalar& factor_ = Scalar(1.0),
                             const Vector& offset_ = Vector::Zero())
        : factor(factor_)
        , offset(offset_)
    {
    }

    void operator()(Vector& vector) const { vector = vector * factor + offset; }

    Scalar factor;
    Vector offset;
  };

  // generator more or less evenly distributed points on a sphere
  template <class Vector>
  void generatePoints(unsigned int N, std::vector<Vector>& points,
                      const Vector& center = Vector::Zero(), typename Vector::Scalar radius = 1.0,
                      unsigned int iterations = 500, typename Vector::Scalar dist = 1e-3,
                      typename Vector::Scalar dt = 0.1)
  {
    typedef typename Vector::Scalar Scalar;
    // construct vector
    std::vector<Vector> newPoints;
    newPoints.reserve(N);

    // generate N points on unit sphere
    std::generate_n(std::back_inserter(newPoints), N, PointOnUnitSphereGenerator<Vector>());
    // optimize these points on unit sphere
    optimizePointsOnUnitSphere(newPoints, iterations, dist, dt);

    // map from unit sphere to radius-sphere at center
    std::for_each(newPoints.begin(), newPoints.end(),
                  ScaleAndTranslateFunctor<Vector>(radius, center));

    // copy the new vectors into the output
    std::copy(newPoints.begin(), newPoints.end(), std::back_inserter(points));
  }
}

#endif // SPHEP_GENERATOR_HH
